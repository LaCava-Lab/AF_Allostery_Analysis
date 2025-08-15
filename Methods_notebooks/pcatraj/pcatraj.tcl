set colorids {0 3 7 1}
set num_frames [molinfo top get numframes]

for {set i 0} {$i < $num_frames} {incr i 1} {
	if {[expr {$i % 5}] == 0} {
    		set color_index [expr {($i/5) % [llength $colorids]}]
    		set colorid [lindex $colorids $color_index]
	}
	mol modcolor 0 top ColorID $colorid
	animate goto $i
    	display update
	# Render the current frame as an image (e.g., with Tachyon)
    	set fname [format "c:/output_frames/frame%04d.tga" $i]
    	render TachyonInternal $fname
}