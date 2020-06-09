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
// separate GFP from DAPI channel
run("Stack to Images");
selectWindow( final_filename+"-0001"); // DAPI
run("Blue");

// select GFP channel
selectWindow( final_filename+"-0002"); // GFP
	
// image processing
run("Subtract Background...", "rolling=8");
run("Set Scale...", "distance=0.7752 known=1 pixel=1 unit=µm global");

// machine learning segementation 
selectWindow( final_filename+"-0002");
run("Trainable Weka Segmentation");
wait(2000);
selectWindow("Trainable Weka Segmentation v3.2.5");
call("trainableSegmentation.Weka_Segmentation.loadClassifier", "select path");
call("trainableSegmentation.Weka_Segmentation.getResult");

// binarize image
selectWindow("Classified image");
run("Threshold...");
setThreshold(1, 255);
run("Convert to Mask");
run("Make Binary");
run("Watershed");
selectWindow("Trainable Weka Segmentation v3.2.5");
close();

// combine DAPI + GFP + segmented image into one stack
run("Images to Stack", "name=Stack title=[] use");
selectWindow("Stack");

saveAs("Tiff", output+"\\"+final_filename +"_stack_with_segmented_cells.tif");
close();

}
setBatchMode(false);