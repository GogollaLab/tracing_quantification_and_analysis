// Choose image, roiset, and where to save it

filepath=File.openDialog("Select a File"); 
output = getDirectory("Choose were to save it");
file=File.openAsString(filepath);
open(filepath);
filename = File.getName(filepath);
section ="";
bregma="";

length = lengthOf(filename)-31;
final_filename = substring(filename,0,length);


// Define Section number and Distance from Bregma
Dialog.create("Define Section ")
Dialog.addString("section number", section);
Dialog.addString("Distance from Bregma", bregma);
Dialog.show();
section = Dialog.getString();
bregma = Dialog.getString();


type= "no";
title ="";
roi ="";
run("Set Measurements...", "area centroid display nan redirect=None decimal=3");


// wait for user, adjust ROIs manually
title_1 = "WaitForUserDemo";
msg = "Adjust ROIs ";
waitForUser(title_1, msg);

run("Stack to Images");
n = nImages; 
list = newArray(n); 
for (i=1; i<=n; i++) { 
       selectImage(i); 
       list[i-1] = getTitle;
       } 

for(k=1; k <=list.length; k++) {
	if (list[k-1] == final_filename+"-0001 largest Hessian eigenvalues"){ 
       	number = "0001";}
    else{
    	number = "0002";
    }
}

selectWindow(final_filename + "-"+number+" largest Hessian eigenvalues");


// find out how many ROIs there are, get the name from each roi; then analyze only as much as there are
n = roiManager("count");
for (i=0; i<n; i++) {
      roiManager("select", i);
      roi_name = call("ij.plugin.frame.RoiManager.getName", i);
      // get roi_area and add to results
      getStatistics(area, mean);
      // run particle analyzer to count cells      
      nBins = 256;
  	  //run("Clear Results");
      //row = 0;
      getHistogram(values, counts, nBins);
      //for (j=0; j<nBins; j++) {
      setResult("Count", i,  counts[0]);
      setResult("Total Area",  i, area);
      setResult("Label",  i, roi_name);
      pixel_density = counts[0]/area;
      setResult("Pixel Density [Pixel/um²]", i , pixel_density);   
        //row++;
      
      }
      updateResults();
   

// save excel file with the data for every single cell
selectWindow("Results");
saveAs("Results", output+"\\" +bregma + "_" +section +"_axon analysis.xls");



selectWindow("Results");
run("Close");
selectWindow("ROI Manager");
run("Close");

run("Close All");





