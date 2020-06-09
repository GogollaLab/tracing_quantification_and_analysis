// Choose image, roiset, and where to save it

filepath=File.openDialog("Select a File"); 
output = getDirectory("Choose were to save it");
file=File.openAsString(filepath);
open(filepath);
filename = File.getName(filepath);
section ="";
bregma="";


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
selectWindow("Classified image");
setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(0, 0);
run("Convert to Mask");


// find out how many ROIs there are, get the name from each roi; analyze only as many as there are
n = roiManager("count");
for (i=0; i<n; i++) {
      roiManager("select", i);
      roi_name = call("ij.plugin.frame.RoiManager.getName", i);
      // get roi_area and add to results
      getStatistics(area, mean);
      // run particle analyzer to count cells      
      run("Analyze Particles...", "size=70-1000 circularity=0.30-1.0 display summarize slice");
      }

// save excel file with the data for every single cell
selectWindow("Results");
saveAs("Results", output+"\\" +bregma + "_" +section + "_" +roi+"_single_cells.xls");

// save summary table; read, fill the results window, save this window
selectWindow("Summary");
text = getInfo("window.contents");
lines = split(text, "\n");
labels = split(lines[0], "\t");
run("Clear Results");
for (k=1; k<lines.length; k++) {
      items=split(lines[k], "\t");
      setResult("Label", k-1, items[0]);
      for (j=1; j<labels.length; j++)
         setResult(labels[j], k-1, items[j]);
}

// change the results table. Replace Total Area value with the true total ROI area and replace cryptic image name with ROI name
n = roiManager("count");
  for (i=0; i<n; i++) {
      roiManager("select", i);
      roi_name = call("ij.plugin.frame.RoiManager.getName", i);
      getStatistics(area, mean);
      setResult("Total Area",  i, area);
      setResult("Label",  i, roi_name);
      }
   updateResults();
selectWindow("Results");
saveAs("Results", output+"\\" +bregma + "_" +section + "_" +roi+"_summary.xls");

selectWindow("Results");
run("Close");
selectWindow("ROI Manager");
run("Close");
selectWindow("Summary");
run("Close");

run("Close All");


	



