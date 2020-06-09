// choose folders for batch process
mainDir = getDirectory("Choose a main directory "); 
mainList = getFileList(mainDir); 
output = getDirectory("Choose were to save it");

setBatchMode(true);
// loop through every file
for (i=0; i<mainList.length; i++) {
	path = mainDir + mainList[i];
	open(path);
	if (bitDepth() != 8){
		run("8-bit");
	}
length = lengthOf(mainList[i])-4;
final_filename = substring(mainList[i],0,length);
run("Remove Slice Labels");

// settings for 1024x1024 confocal images
run("Set Scale...", "distance=1 known=1.4 pixel=1 unit=µm global");

// settings for 512x512 confocal images 
//run("Set Scale...", "distance=0.3303 known=1 pixel=1 unit=µm");

// separate GFP from DAPI channel
run("Stack to Images");
selectWindow( final_filename+"-0002");
run("Enhance Contrast...", "saturated=1.0");


// select GFP channel
selectWindow( final_filename+"-0001");

// Hessian ridge detection to segment axons out
run("FeatureJ Hessian", "largest absolute smoothing=0.5");

//thresholding image Isodata 0,5
selectWindow( final_filename+"-0001 largest Hessian eigenvalues");
run("8-bit");
setAutoThreshold("IsoData");
run("Threshold...");
setThreshold(0, 25);
run("Convert to Mask");
run("Make Binary");

// increase contrast in original gfp image
selectWindow( final_filename+"-0001");
run("Enhance Contrast...", "saturated=0.3");


// combine DAPI + GFP + segmented image into one stack
run("Images to Stack", "name=Stack title=[] use");
selectWindow("Stack");

saveAs("Tiff", output+"\\"+final_filename +"_stack_with_segmented_axons.tif");
close();


}
setBatchMode(false);