"""compare_regression_rds.py: independent reader for the two .rds files that
check_dgt_060_regression.R writes (lists of data.frames of posterior draws),
parsing R's XDR serialization directly (versions 2 and 3, gzip or plain) and
applying the script's own criterion, max |diff| < 1e-10 per column. Written so
the comparison can be rerun, and its record regenerated, on a machine without
R. Usage:
  python3 -I verification/compare_regression_rds.py \
      verification/dgt_regression_installed.rds verification/dgt_regression_devel.rds \
      [commit note]
"""

import gzip, struct, sys, statistics

class R:
    def __init__(self, b):
        self.b = b; self.p = 0; self.refs = []
    def i(self):
        v = struct.unpack(">i", self.b[self.p:self.p+4])[0]; self.p += 4; return v
    def d(self, n):
        v = struct.unpack(">%dd" % n, self.b[self.p:self.p+8*n]); self.p += 8*n; return list(v)
    def length(self):
        n = self.i()
        if n == -1:
            hi = self.i(); lo = self.i(); n = (hi << 32) + lo
        return n
    def item(self):
        flags = self.i()
        t = flags & 0xFF; has_attr = flags & 0x200; has_tag = flags & 0x400
        if t == 254: return None                      # NILVALUE
        if t == 255: return self.refs[(flags >> 8) - 1] # REFSXP
        if t == 242: return None                      # EMPTYENV etc. not expected
        if t == 1:                                    # SYMSXP
            s = self.item(); self.refs.append(s); return s
        if t == 9:                                    # CHARSXP
            n = self.i()
            if n == -1: return None
            s = self.b[self.p:self.p+n].decode("utf-8", "replace"); self.p += n; return s
        if t in (2, 6, 17):                            # LISTSXP / LANGSXP / DOTSXP
            attr = self.item() if has_attr else None
            tag = self.item() if has_tag else None
            car = self.item(); cdr = self.item()
            out = [(tag, car)]
            if isinstance(cdr, list): out += cdr
            return out
        if t == 10:                                   # LGLSXP
            n = self.length(); v = [self.i() for _ in range(n)]
        elif t == 13:                                 # INTSXP
            n = self.length(); v = [self.i() for _ in range(n)]
        elif t == 14:                                 # REALSXP
            n = self.length(); v = self.d(n)
        elif t == 16:                                 # STRSXP
            n = self.length(); v = [self.item() for _ in range(n)]
        elif t == 19:                                 # VECSXP
            n = self.length(); v = [self.item() for _ in range(n)]
        else:
            raise ValueError("unhandled SEXP type %d at %d" % (t, self.p))
        attrs = self.item() if has_attr else None
        names = None
        if attrs:
            for tag, val in attrs:
                if tag == "names": names = val
        if names is not None and t == 19:
            return dict(zip(names, v))
        if names is not None:
            return dict(zip(names, v))
        return v

def read_rds(path):
    raw = open(path, "rb").read()
    if raw[:2] == b"\x1f\x8b": raw = gzip.decompress(raw)
    assert raw[:2] == b"X\n", "not XDR"
    r = R(raw); r.p = 2
    ver = r.i(); r.i(); r.i()
    if ver == 3:
        n = r.i(); r.p += n
    return r.item()

a = read_rds(sys.argv[1]); b = read_rds(sys.argv[2])
note = sys.argv[3] if len(sys.argv) > 3 else ""
print("Comparison of %s (dgt %s) and %s (dgt %s)%s, produced by" % (
    sys.argv[1], a["version"][0], sys.argv[2], b["version"][0], (", " + note) if note else ""))
print("verification/compare_regression_rds.py with the criterion of")
print("check_dgt_060_regression.R, max |diff| < 1e-10 per column.")
print("Medians of the saved draws (second file):\n")
for nm in [k for k in a if k != "version"]:
    d = b[nm]
    print("  %-18s median ICC_Y_abs %.4f  rel %.4f  eta %.4f  I %.4f" % (nm,
        statistics.median(d["icc_Y_abs"]), statistics.median(d["icc_Y_rel"]),
        statistics.median(d["icc_eta"]), statistics.median(d["icc_I"])))
print("\ninstalled %s vs devel %s" % (a["version"][0], b["version"][0]))
ok = True
for nm in [k for k in a if k != "version"]:
    for col in a[nm]:
        x, y = a[nm][col], b[nm][col]
        if len(x) != len(y):
            print("%-18s %-10s LENGTH DIFFERS" % (nm, col)); ok = False; continue
        dmax = max(abs(p - q) for p, q in zip(x, y))
        flag = "identical" if dmax < 1e-10 else "DIFFERS (max |diff| = %.3g)" % dmax
        if dmax >= 1e-10: ok = False
        print("%-18s %-10s %s" % (nm, col, flag))
print("\nAll draws identical between the two versions." if ok else "\nDIFFERENCES FOUND.")
