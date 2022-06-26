 set term postscript color enhanced
 set size 0.9,1
 rearth=   6371.0000000000000     
 au=   149597870.69999999     
 rpau=rearth/au
 set title "Itokawa"
 set key bottom right
 set output "rearth_sigma.eps"
 set arrow from  -3.0000000000000000      ,1 to    3.0000000000000000      ,1 nohead
 set xlabel "{/Symbol s}_{LOV}"
 set ylabel "Min. distance [R_{Earth}]"
 p[  -3.0000000000000000      :   3.0000000000000000      ][:] "Radau_sml.dmi"  u (  -6.0000000000000000      *(          54 -$1)/         100 ):($8/rpau)tit "hyperbolic" w lp ps 0.9, "" u (  -6.0000000000000000      *(          54 -$1)/         100 ):($5/rpau) tit "cubic spline" w lp ps 0.5
 set output "dum.eps"
 set term x11
