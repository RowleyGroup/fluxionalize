#! /usr/bin/tclsh

set psf "mol_wb.psf"
set outfile [open mol_cluster.dat w]
set outfil2 [open mol_rmsdtt.dat w]
set w1 7
set w2 30
set title [format "%-*s" $w1 "frame"]
for {set j 0} {$j <= 3} {incr j} {
 set dcd "output/$j/mol.job0.$j.dcd"
 mol load psf $psf dcd $dcd
 set nf [molinfo top get numframes]
 set frame0 [atomselect top "not water" frame 0]
 set sel [atomselect top "not water"]
 append title [format "%-*s" $w2 "mol$j"]
#measures RMSD for each frame aligning to frame 0
 for {set i 0} {$i <= $nf} {incr i} {
  $sel frame $i
  $sel move [measure fit $sel $frame0]
  append outline($i) [format "%-*s" $w2 "[measure rmsd $sel $frame0]"]
 }
#finds clusters for each replica, writing a temp dcd for each frame and then combining them to one dcd per cluster
 set clust [measure cluster $sel num 5 distfunc rmsd]
 puts $outfile $clust
 set h 0
 foreach cluster $clust {
   set outclust [open clust$j.$h.txt w]
   foreach frame $cluster {
    puts $outclust $frame
    }
  }
  exec ./catdcd -o clusters/cluster$j.$h.dcd  -f $outclust $dcd
  incr h
  close $outclust
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

