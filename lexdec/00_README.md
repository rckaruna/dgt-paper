# Illustration 1 (Section 9.1): lexical decision reaction times

`01_lexdec.R` fits the lognormal crossed random-effects model of Section 9.1
to the correct-response reaction times of `languageR::lexdec` (21 subjects,
79 words) with `brms` and computes the link-scale and response-scale ICCs,
the overestimation ratio, and the D-study with `dgt`. The numbers quoted in
Section 9.1 are printed under `HEADLINE` at the end of `log_01.txt`.

    Rscript 01_lexdec.R 2>&1 | tee log_01.txt

The illustration was part of the submitted version of the paper; this
script was added in the revision so that the illustration is reproducible
like the others. The fit is cached as `lexdec_fit.rds` (not tracked).
