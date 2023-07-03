# Contact analysis
#
# prints output to 1 files:
# contact_all.dat
# 1 line per contact and frame
# frame number, contacting residues, minimal distance
# example:
# protein residues 1 (ARG), 5 (LYS) are in contact with 10 (CL)
# 1 1 ARG - 10 CL 3.50000
# 1 5 LYS - 10 CL 2.90000

source bigdcd.tcl
proc fResPair { frame } {
  global sel1 sel2 all 
  # resid map for every atom
  set allResid [$all get residue]
  # resname map for every atom
  set allResname [$all get resname]
  # create resid->resname map
  foreach resID $allResid resNAME $allResname {
    set mapResidResname($resID) $resNAME
  }	
	
  # output file for printing out all contacts
  set fall [open "contact_all.dat" a+]
  
  # position for every atom
  set allPos [$all get {x y z}]
  set cutoff cutoff_dist 
	 
  # extract the pairs. listA and listB hold corresponding pairs from selections 1 and 2, respectively.
  # To run foreach loop over a list for variable assignment and not need any 'body': use 'break' 
  # measure Contacts: Returns two lists of atom indices
  foreach {listA listB} [measure contacts $cutoff $sel1 $sel2] break
  
 # go through the pairs, assign distances
    foreach indA $listA indB $listB {
      # calculated distance between 2 atoms
      set dist [vecdist [lindex $allPos $indA] [lindex $allPos $indB]]

      # get information about residue id's
      set residA [lindex $allResid $indA]
      set residB [lindex $allResid $indB]
      # following code can be uncommented for testing purposes
      #set resnameA [lindex $allResname $indA]
      #set resnameB [lindex $allResname $indB]
      #puts "$residA $resnameA $residB $resnameB $dist"
    
      # results will be stored in [residA,residB][list of distances] array
      lappend contactTable($residA,$residB) $dist
    }
	
	foreach name [array names contactTable] {
      # assign residue names to residue numbers
      foreach {residA residB} [split $name ,] break
      foreach {tmp resnameA} [split [array get mapResidResname $residA] ] break
      foreach {tmp resnameB} [split [array get mapResidResname $residB] ] break
      # get minimum contact distance for the pair
      foreach {tmp distanceList} [array get contactTable $name] break
      set minDistance [lindex [lsort -real $distanceList] 0]

      # print to output file
      puts $fall [format "%-5d %5d %3s - %5d %3s - %f" $frame $residA $resnameA $residB $resnameB $minDistance]
      flush $fall
	 
  }
   # delete the contact table - will be created again in the next loop
    if {[info exists contactTable]} {
      unset contactTable
    }
  close $fall
  
 }
  
  set mol [mol new top_file waitfor all]
  set all [atomselect top all]
  
  # create specified atom selections
  set sel1 [atomselect top "selection1"]
  set sel2 [atomselect top "selection2"]
  
  #load dcd
  bigdcd fResPair dcd_file
  bigdcd_wait
  
quit
