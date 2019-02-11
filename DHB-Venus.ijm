macro "DHB-Ven" {
	firstLoop = true;
	
	//Identification of folders
	//create an array containing the names of the files in the directory path
	dir=getDirectory("Choose a Directory");
	dirParent=File.getParent(dir);
	dirName=File.getName(dir);
	outputPath=dirParent+File.separator+dirName+"_OUT";
	list=getFileList(dir);
	Array.sort(list);
	numberOfFolders=0;
	File.makeDirectory(outputPath);
	
	//count the number of folders
	for (i=0; i<list.length; i++) {
		if (endsWith(list[i], ".oif")==true) {		
			numberOfFolders++;
		}
	}
	
	//check that the directory contains folders
	if (numberOfFolders==0) {
		beep();
		exit("No folders")
	}
	
	//create a an array containing only the names of the folders in the directory path
	folderArray=newArray(numberOfFolders);
	count=0;
	for (i=0; i<list.length; i++) {
		if (endsWith(list[i], ".oif")) {
			folderArray[count]=list[i];
			count++;
		}
	}

	//set parameters
	enhanceContrastOptions=newArray("0", "0.1", "0.2", "0.3", "0.4", "None");
	thresholdMethods=getList("threshold.methods");
	Dialog.create("Segmentation parameters");
	Dialog.setInsets(0, 115, 0);
	Dialog.addMessage("NUCLEI");
	Dialog.addNumber("Rolling ball", 30);
	Dialog.addToSameRow();
	Dialog.addChoice("Enhance Contrast", enhanceContrastOptions, enhanceContrastOptions[2]);
	Dialog.addNumber("Mean filter", 1);
	Dialog.addToSameRow();
	Dialog.addChoice("setAutoThreshold", thresholdMethods, "Triangle");
	Dialog.addSlider("Watershed solidity", 0, 1, 0.89);
	Dialog.addNumber("Open", 4);
	Dialog.addToSameRow();
	Dialog.addNumber("Erode", 1);
	Dialog.addNumber("Size (min)", 20);
	Dialog.addToSameRow();
	Dialog.addNumber("Size (max)", 200);
	Dialog.show();
	rollingBall=Dialog.getNumber();
	normalize=Dialog.getChoice();
	mean=Dialog.getNumber();
	threshold=Dialog.getChoice();
	wSol=Dialog.getNumber();
	openI=Dialog.getNumber();
	erodeI=Dialog.getNumber();
	min=Dialog.getNumber();
	max=Dialog.getNumber();


	for (i=0; i<folderArray.length; i++) {
		run("Bio-Formats", "open=["+dir+"\\"+folderArray[i]+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		run("Stack to Images");
		for (j=1; j<=nImages; j++) {
			selectImage(j);
			run("Grays");
			title=getTitle();
			for (k=1; k<=nImages; k++) {
				channel=d2s(k, 0);
				if(indexOf(title, "c:"+channel)!=-1) {
					rename("s_C00"+j);
				}
			}
		}

		dapiSegmentation(rollingBall, normalize, mean, threshold, wSol, openI, erodeI, min, max);
		venusSegmentation();
		venusAndDapi();
		voronoi();
		
		//Nucleus by nucleus analysis
		selectImage("DHBVen-AND-DAPI");
		run("Analyze Particles...", " size=0-Infinity show=[Count Masks] display clear");
		rename("DAPI-CountMasks");
		run("8-bit");
		close("DHBVen-AND-DAPI");
		close("s_C001-EEA");
		localNucleiMean=newArray(nResults);
		localNucleiIntDen=newArray(nResults);
		localNucleiRawIntDen=newArray(nResults);
		localNucleiArea=newArray(nResults);
		localCytoplasmMean=newArray(nResults);
		localCytoplasmIntDen=newArray(nResults);
		localCytoplasmRawIntDen=newArray(nResults);
		localCytoplasmArea=newArray(nResults);
		localDataname=newArray(nResults);
		masksNumber=nResults;
		for (m=1; m<=masksNumber; m++) {
			nameNucleous="DAPI-"+m;
			selectImage("DAPI-CountMasks");
			run("Duplicate...", "title="+nameNucleous);
			setThreshold(m, m);
			run("Convert to Mask");
			run("Make Binary");
			run("BinaryReconstruct ", "mask=CytoplasmMasks seed="+nameNucleous+" create white");
			rename("CellMask"+"-"+m);
			imageCalculator("XOR create", nameNucleous,"CellMask"+"-"+m);
			rename("CytoplasmMask"+"-"+m);
			close("CellMask"+"-"+m);

			//Nuclei: analyze mean gray value, integrated density and area
			run("Set Measurements...", "area mean integrated redirect=s_C002 decimal=2");
			selectImage(nameNucleous);
			run("Analyze Particles...", "display exclude clear");
			localNucleiMean[m-1]=getResult("Mean", 0);
			localNucleiIntDen[m-1]=getResult("IntDen", 0);
			localNucleiRawIntDen[m-1]=getResult("RawIntDen", 0);
			localNucleiArea[m-1]=getResult("Area", 0);
			localDataname[m-1]=folderArray[i];
			close(nameNucleous);

			//Cytoplasms: analyze mean gray value, integrated density and area
			selectImage("CytoplasmMask"+"-"+m);
			run("Analyze Particles...", "display exclude clear");
			
			//Selection of the largest particle as cytoplasm mask
			iteLargerArea=0;
			if (nResults > 1) {
				localCytoplasmArea[m-1]=getResult("Area", 0);
				for (n=1; n<nResults; n++) {
					areaFraction=getResult("Area", n);
					if (areaFraction > localCytoplasmArea[m-1]) {
						localCytoplasmArea[m-1] = areaFraction;
						iterLargerArea = n;
					}
				}
			} else if (nResults < 1) {
				localCytoplasmArea[m-1]=-1;
			} else {
				localCytoplasmArea[m-1]=getResult("Area", 0);
			}

			//the cytoplasm is taken into account only if its area is equal to or greater than the nucleus area
			if (localCytoplasmArea[m-1] >= localNucleiArea[m-1]) {
				localCytoplasmMean[m-1]=getResult("Mean", iteLargerArea);
				localCytoplasmIntDen[m-1]=getResult("IntDen", iteLargerArea);
				localCytoplasmRawIntDen[m-1]=getResult("RawIntDen", iteLargerArea);
			} else {
				localCytoplasmMean[m-1]=-1;
				localCytoplasmIntDen[m-1]=-1;
				localCytoplasmRawIntDen[m-1]=-1;
				localCytoplasmArea[m-1]=-1;
			}
			//Close "CytoplasmMask"+"-"+m
			close("CytoplasmMask"+"-"+m);
		}
		
		//Store results
		if (firstLoop) {
			nucleiMean=localNucleiMean;
			nucleiIntDen=localNucleiIntDen;
			nucleiRawIntDen=localNucleiRawIntDen;
			nucleiArea=localNucleiArea;
			cytoplasmMean=localCytoplasmMean;
			cytoplasmIntDen=localCytoplasmIntDen;
			cytoplasmRawIntDen=localCytoplasmRawIntDen;
			cytoplasmArea=localCytoplasmArea;
			dataname=localDataname;
			firstLoop=false;
		} else {
			nucleiMean=Array.concat(nucleiMean, localNucleiMean);
			nucleiIntDen=Array.concat(nucleiIntDen, localNucleiIntDen);
			nucleiRawIntDen=Array.concat(nucleiRawIntDen, localNucleiRawIntDen);
			nucleiArea=Array.concat(nucleiArea, localNucleiArea);
			cytoplasmMean=Array.concat(cytoplasmMean, localCytoplasmMean);
			cytoplasmIntDen=Array.concat(cytoplasmIntDen, localCytoplasmIntDen);
			cytoplasmRawIntDen=Array.concat(cytoplasmRawIntDen, localCytoplasmRawIntDen);
			cytoplasmArea=Array.concat(cytoplasmArea, localCytoplasmArea);
			dataname=Array.concat(dataname, localDataname);
		}
		
		//Clean up
		cleanUp();
	}

	//Results table
	title1 = "Results table";
	title2 = "["+title1+"]";
	f = title2;
	run("Table...", "name="+title2+" width=500 height=500");
	print(f, "\\Headings:n\tdataname\tNucleus Mean\tNucleus IntDen\tNucleus RawIntDen\tNucleus RawIntDen / Area\tNucleus Area\tCytoplasm Mean\tCytoplasm IntDen\tCytoplasm RawIntDen\tCytoplasm RawIntDen / Area\tCytoplasm Area");
	for (i=0; i<nucleiMean.length; i++) {
		n=i+1;
		ridArea=nucleiRawIntDen[i]/nucleiArea[i];
		criArea=cytoplasmRawIntDen[i]/cytoplasmArea[i];
		print(f, n + "\t" + dataname[i] + "\t" + nucleiMean[i] + "\t" + nucleiIntDen[i] + "\t" + nucleiRawIntDen[i] + "\t" + ridArea + "\t" + nucleiArea[i] + "\t" + cytoplasmMean[i] + "\t" + cytoplasmIntDen[i] + "\t" + cytoplasmRawIntDen[i] + "\t" + criArea + "\t" + cytoplasmArea[i]);
	}

	saveAs("txt", outputPath+File.separator+"Results table - DHB-Venus");
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Functions

	function dapiSegmentation (a, b, c, d, e, f, g, h, i) {
		selectImage("s_C001");
		run("8-bit");
		run("Subtract Background...", "rolling="+a);
		run("Enhance Contrast...", "saturated="+b+" normalize");
		run("Mean...", "radius="+c);
		setAutoThreshold(d+" dark");
		run("Convert to Mask");
		run("Set Measurements...", "  redirect=None decimal=2");
		run("Analyze Particles...", "  show=Masks display exclude clear");
		rename("s_C001-EE");
		watershedSolidity(e);
		run("Options...", "iterations=1 count="+f+" do=Open");
		run("Options...", "iterations=1 count="+g+" do=Erode");
		run("Analyze Particles...", "size="+h+"-"+i+" show=Masks display exclude clear");
		rename("s_C001-EEA");
		close("s_C001");
		close("Result");
	}

	function venusSegmentation () {
		selectImage("s_C002");
		run("Duplicate...", "title=venus-segmentation");
		run("8-bit");
		setAutoThreshold("Triangle dark");
		run("Convert to Mask");
		run("Options...", "iterations=3 count=1 do=Close");
		run("Fill Holes");
		run("Set Measurements...", "  redirect=None decimal=2");
		run("Analyze Particles...", "size=15-Infinity show=Masks display clear");
		rename("Venus-Masks");
		close("venus-segmentation");
	}

	function venusAndDapi () {
		imageCalculator("AND create", "s_C001-EEA","Venus-Masks");
		run("Fill Holes");
		rename("DHBVen-AND-DAPI");
	}

	function voronoi () {
		run("Duplicate...", "title=territories");
		run("Voronoi");
		setThreshold(1, 255);
		run("Convert to Mask");
		run("Invert");
		imageCalculator("AND create", "territories","Venus-Masks");
		rename("CytoplasmMasks");
		run("Make Binary");
		close("territories");
		close("Venus-Masks");
	}
	
	function watershedSolidity (maxSolidity) {
		rename("Original");
		run("Set Measurements...", "shape stack redirect=None decimal=2");
		run("Analyze Particles...", "display exclude clear record");
		run("Classify Particles", "class[1]=Solidity operator[1]=<= value[1]="+maxSolidity+
		" class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]="+
		"-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
		rename("Watershed");
		imageCalculator("XOR create", "Original","Watershed");
		rename("No-watershed");
		selectImage("Watershed");
		run("Watershed");
		imageCalculator("OR create", "Watershed","No-watershed");
		rename("Result");
		close("Original");
		close("Watershed");
		close("No-watershed");
		selectWindow("Result");
	}
	
	function cleanUp() {
		requires("1.30e");
		if (isOpen("Results")) {
			selectWindow("Results");
			run("Close");
		}
		if (isOpen("Threshold")) {
			selectWindow("Threshold"); 
			run("Close");
		}
		if (isOpen("ROI Manager")) {
			selectWindow("ROI Manager");
			run("Close");
		}   
		if (isOpen("Channels")) {
			selectWindow("Channels");
			run("Close");
		}
		if (isOpen("Log")) {
			selectWindow("Log");
			run("Close");
		}
		while (nImages()>0) {
			selectImage(nImages());  
			run("Close");
		}
	}
}