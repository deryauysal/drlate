



logit e401 incst age25st age25sqst marr fsize
predict nix
gen w1 = e401/nix
gen w0 = (1-e401)/(1-nix)

ge y1hat =  net_tfa*w1
ge y0hat =  net_tfa*w0

ge numhat=y1hat-y0hat


ge d1hat =  p401*w1
ge d0hat =  p401*w0

ge denomhat=d1hat-d0hat


sum numhat denomhat
