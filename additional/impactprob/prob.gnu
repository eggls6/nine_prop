 set term postscript color enhanced
 set size 0.9,1
 rearth=   6371.0000000000000     
 au=   149597870.69999999     
 rpau=rearth/au
 set title "Itokawa"
 set key bottom right
 set output "d_prob.eps"
 set xlabel "Probability"
 set xlabel "Minimum distance [R_{Earth}]"
 p[:][0:1] "prob.dat"  u ($2/rpau):3 notit w l lw 3
 set output "dum.eps"
 set term x11
