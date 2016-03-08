#! /usr/bin/tclsh

set psf "mol_wb.psf"
##Must set catdcd location manually, see example below##
set catdcd "/home/rowley_group/lib/vmd/plugins/LINUXAMD64/bin/catdcd5.1/catdcd"
set outfile [open mol_cluster.dat w]
set outfil2 [open mol_rmsdtt.dat w]
set w1 7
set w2 30
set title [format "%-*s" $w1 "frame"]
for {set j 0} {$j <= 3} {incr j} {
 set dcd "output/$j/mol.job0.$j.dcd"
 mol load psf $psf dcd $dcd
 set nf [molinfo top get numframes]
 set frame0 [atomselect top "not water" frame 1000]
 set sel [atomselect top "not water"]
 append title [format "%-*s" $w2 "mol$j"]
#measures RMSD for each frame aligning to frame 1000
 for {set i 1000} {$i <= $nf} {incr i} {
  $sel frame $i
  $sel move [measure fit $sel $frame0]
  append outline($i) [format "%-*s" $w2 "[measure rmsd $sel $frame0]"]
 }
#finds clusters for each replica, writing a temp dcd for each frame and then combining them to one dcd per cluster
 set clust [measure cluster $sel num 5 distfunc rmsd]
 puts $outfile $clust
 set h 0
 foreach cluster $clust {
   set k 0
   foreach frame $cluster {
    animate write dcd temp$k.dcd beg $frame end $frame sel $sel waitfor all top
    incr k
    }
   set nf [llength $cluster]
   for {set l 0} {$l < [expr $nf-1]} {incr l} {
    set m [expr $l+1]
    exec $catdcd -o temp.dcd temp$l.dcd temp$m.dcd
    exec mv temp.dcd temp$m.dcd
    exec rm temp$l.dcd
  }
  exec mv temp$m.dcd clusters/cluster$j.$h.dcd
  incr h
 }
}
puts $outfil2 $title
for {set h 1} {$h <=$nf} {incr h} {
 set newline [format "%-*s" $w1 $h]
 puts $outfil2 "$newline $outline($h)"
}
close $outfile
close $outfil2
exit
EOF
