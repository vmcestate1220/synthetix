load AF-Q43127-F1.pdb, gs2

select cat_domain, resi 111-430

set bg_rgb, [0.0863, 0.1059, 0.1333]
hide everything
show cartoon, cat_domain
set cartoon_transparency, 0.0
set ray_shadows, 0
set ambient, 0.3
set specular, 0.15
set ray_trace_mode, 1
set ray_trace_color, black
set ray_opaque_background, 1
set antialias, 2

color 0xFF7D45, gs2
color 0xFFDB13, gs2 and b > 50
color 0x65CBF3, gs2 and b > 70
color 0x0053D6, gs2 and b > 90

orient cat_domain
zoom cat_domain, 3

ray 1600, 1000
png GS2_AlphaFold_Structure.png, dpi=200
