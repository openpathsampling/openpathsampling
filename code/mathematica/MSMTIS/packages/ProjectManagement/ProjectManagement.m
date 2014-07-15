(* ::Package:: *)

(* ProjectManagement Package *)
(* Version 0.1 *)

(* Last edited on Jan 24th, 2010 *)

BeginPackage["ProjectManagement`"];


(* Messages *)

InitializeProject::usage = "Creates aproject folder if necessary and sets all path variables correctly.";

$projectName::usage = "Name of the project";
$projectPath::usage = "Path to the project folder";
$imagePath::usage = "Path to folder for all images";
$resourcePath::usage = "Path to ressource folder";
$literaturePath::usage = "Path to literature folder includeing bibTeX files.";
$manuscriptPath::usage = "Path to manuscript folder";
$backupPath::usage = "Path to backup folder";
$localPackagesPath::usage = "Path to package folder";

$DefaultImageResolution::usage = "Standard Image Resolution for exporting graphics";

ToImageFile::usage = "Add the correct image path to a string file given.";
ToResourceFile::usage = "Add the correct resource path to a string file given.";
ToLiteratureFile::usage = "Add the correct literature path to a string file given.";
ToManuscriptFile::usage = "Add the correct manuscript path to a string file given.";

CopyPackagesToLocalFolder::usage = "Copies all non kernel packages to a local folder packages";
LoadLocalPackages::usage = "Loads all packages in the local packages folder";

ProjectSaveBackUp::usage = "Stores a copy of the current file in the backup folder with the current date attached to the filename.";

ListResources::usage = "ListResources[expression] returns a list of all resources that match expression";
SaveGraphics::usage = "SaveGraphics[names,graphics] saves all graphics with their respective names in the graphics path.";
SaveGraphics::res = "Using resolution of `1` dpi.";
CloneProject::usage = "Copy the current project to another location and name.";
ExportProject::usage = "Export parts of the folders to a zip file in the folder exports. Without options Folders given the images and manuscipt folder will be exported.";

CreateSaveGraphicsPalette::usage = "Opens a palette to conveniently save all graphics to the image folder";
SaveGraphicsByName::usage = "Save a graphics object with the name given.";
SaveAllGraphicsByName::usage = "Saves all graphics object starting with the given prefix";
OpenProjectFolder::usage = "Opens the project folder in the Mac Finder.";
OpenImageFolder::usage = "Open the image folder in the Mac Finder.";

GetFromAddressBook::usage = "Reads an address from the AddressBook";
CreateAuthor::usage = "Creates a userblock from the data in the AddressBook";
GraphicsType::usage = "Option for SaveGraphics[] to specify the graphics types to be exported as a list of strings used in Export[]";

Resolution::usage = "";
Printer::usage = "";
$DefaultImageResolution = 150;
$DefaultGraphicsType = {"pdf","png"};
$saveGraphicsPaletteWindow::usage = "SaveGraphics Palette Object";
$SaveGraphicsPrefix="new_";
$SaveGraphicsContext="gPlot";

Begin["`Private`"];

(* Functions to handle projects *)

InitializeProject[]:=Module[{packagesToBeCopied,$oldPackagesPath,$oldFileName,$basePath,$fileName,$projectNameSelected,$newProjectName},
	If[FileNameTake[NotebookFileName[],{-2,-2}]==FileBaseName[NotebookFileName[]],
		(* Folder has same name as file, that means it has a project folder 
		 * ==> Just get all project Variables and do nothing else !
		 *)
		$projectNameSelected=NotebookFileName[];
		$projectName=FileBaseName[$projectNameSelected];
		$basePath=FileNameTake[$projectNameSelected,{1,-3}];
		$projectPath=FileNameJoin[{$basePath,$projectName}];
		$fileName=FileNameJoin[{$projectPath,$projectName<>".nb"}];
		$imagePath=FileNameJoin[{$projectPath,"images"}];
		$resourcePath=FileNameJoin[{$projectPath,"resources"}];
		$literaturePath=FileNameJoin[{$projectPath,"literature"}];
		$manuscriptPath=FileNameJoin[{$projectPath,"manuscript"}];
		$backupPath=FileNameJoin[{$projectPath,"backup"}];
		$localPackagesPath = FileNameJoin[{$projectPath,"packages"}];
	,
		(* No Project folder seems to exist, ask for a path to create the folder (in form of a .project file)
		 * and copy the file there, rename it and create all useful folder and set the project path variables 
		 *)
		$newProjectName=If[NotebookFileName[]=!=$Failed,FileBaseName[NotebookFileName[]],"NewProject"];
		$projectNameSelected=SystemDialogInput["FileSave",{$newProjectName<>".project",{".project"->{".project"}}},WindowTitle->"Create New Project"];
		If[$projectNameSelected=!=$Canceled,
			$oldFileName=NotebookFileName[];
			$projectName=FileBaseName[$projectNameSelected];
			$basePath=FileNameTake[$projectNameSelected,{1,-2}];
			$projectPath=FileNameJoin[{$basePath,$projectName}];
			$fileName=FileNameJoin[{$projectPath,$projectName<>".nb"}];
			$imagePath=FileNameJoin[{$projectPath,"images"}];
			$resourcePath=FileNameJoin[{$projectPath,"resources"}];
			$literaturePath=FileNameJoin[{$projectPath,"literature"}];
			$manuscriptPath=FileNameJoin[{$projectPath,"manuscript"}];
			$backupPath=FileNameJoin[{$projectPath,"backup"}];
			$localPackagesPath = FileNameJoin[{$projectPath,"packages"}];

			CreateDirectory[$projectPath];
			CreateDirectory[$imagePath];
			CreateDirectory[$resourcePath];
			CreateDirectory[$literaturePath];
			CreateDirectory[$manuscriptPath];
			CreateDirectory[$backupPath];

			CreateDirectory[$localPackagesPath];

			(* If there were already packages present copy these*)

			$oldPackagesPath = FileNameJoin[{FileNameTake[$oldFileName,{1,-2}],"packages"}];
			If[FileExistsQ[$oldPackagesPath],
				packagesToBeCopied = FileNames["*.m*",$oldPackagesPath];
				Map[CopyFile[#,FileNameJoin[{$localPackagesPath,FileNameTake[#]}]]&, packagesToBeCopied];
			];			

			NotebookSave[SelectedNotebook[],$fileName];
			DeleteFile[$oldFileName];
		]
	];
]

CloneProject[]:=Module[{$oldFileName,$basePath,$fileName,$projectNameSelected,$newProjectName},
	If[FileNameTake[NotebookFileName[],{-2,-2}]==FileBaseName[NotebookFileName[]],
		(* Folder has same name as file, that means it has a project folder 
		 * ==> Just get all project Variables and do nothing else !
		 *)
		$projectNameSelectedOld = NotebookFileName[];
		$basePathOld          = FileNameTake[$projectNameSelectedOld,{1,-3}];
		$projectNameOld       = FileBaseName[$projectNameSelectedOld];
		$projectPathOld       = FileNameJoin[{$basePathOld,$projectNameOld}];
		$fileNameOld          = FileNameJoin[{$projectPathOld,$projectNameOld<>".nb"}];

		(* No Project folder seems to exist, ask for a path to create the folder (in form of a .project file)
		 * and copy the file there, rename it and create all useful folder and set the project path variables 
		 *)
		$newProjectName=If[NotebookFileName[]=!=$Failed,FileBaseName[NotebookFileName[]],"NewProject"];
		$projectNameSelected=SystemDialogInput["FileSave",{$newProjectName<>".project",{".project"->{".project"}}},WindowTitle->"Clone Project"];
		If[$projectNameSelected=!=$Canceled,
			$oldFileName=NotebookFileName[];
			$projectName=FileBaseName[$projectNameSelected];
			$basePath=FileNameTake[$projectNameSelected,{1,-2}];
			$projectPath=FileNameJoin[{$basePath,$projectName}];
			$fileName=FileNameJoin[{$projectPath,$projectName<>".nb"}];
			$imagePath=FileNameJoin[{$projectPath,"images"}];
			$resourcePath=FileNameJoin[{$projectPath,"resources"}];
			$literaturePath=FileNameJoin[{$projectPath,"literature"}];
			$manuscriptPath=FileNameJoin[{$projectPath,"manuscript"}];
			$backupPath=FileNameJoin[{$projectPath,"backup"}];
			$localPackagesPath = FileNameJoin[{$projectPath,"packages"}];

			(* Copy folders *)

			CopyDirectory[$projectPathOld, $projectPath];

			(* Move notebook with old name to backup folder *)
			dateStamp=DateString[{"YearShort","Month","Day","_","Hour","Minute","Second"}];
			$backupName=FileNameJoin[{$backupPath,$projectNameOld<>"_"<>dateStamp<>".nb"}];

			RenameFile[FileNameJoin[{$projectPath,$projectNameOld<>".nb"}],$backupName];

			NotebookSave[SelectedNotebook[],$fileName];
		]
	];
]

ToImageFile[name_] := FileNameJoin[{$imagePath,name}]
ToResourceFile[name_] := FileNameJoin[{$resourcePath,name}]
ToLiteratureFile[name_] := FileNameJoin[{$literaturePath,name}]
ToManuscriptFile[name_] := FileNameJoin[{$manuscriptPath,name}]

(* Copy Packages to *)
CopyPackagesToLocalFolder[] := Module[{listOfLoadedPackages,listOfPackagesToBeCopied,result,existingFiles},
	listOfLoadedPackages = SystemInformation["Kernel","InitializationFiles"];
	$localPackagesPath = FileNameJoin[{$projectPath,"packages"}];
	$packagePathsToBeExcluded={"Kernel","Parallel","Paclets"};
	listOfPackagesToBeCopied=Select[listOfLoadedPackages,StringFreeQ[DirectoryName[#],$packagePathsToBeExcluded]&];
	CreateDirectory[$localPackagesPath];
	result=Map[
			CopyFile[#,FileNameJoin[{$localPackagesPath,FileNameTake[#]}]]&	
		,
			listOfPackagesToBeCopied
		];
	existingFiles=Map[FileExistsQ[FileNameJoin[{$localPackagesPath,FileNameTake[#]}]]&,listOfPackagesToBeCopied];
	Apply[And,existingFiles]
];

(* Load All Packages in the local packages folder *)

LoadLocalPackages[] := Module[{packagesToBeLoaded,result},
	If[FileExistsQ[$localPackagesPath],
		packagesToBeLoaded = FileNames["*.m*",$localPackagesPath];
		result=Map[Get, packagesToBeLoaded];
		Transpose[{packagesToBeLoaded,result}]
	];
]

(* BackUp Functionality *)

ProjectSaveBackUp[] := Module[{$fileName, dateStamp, $backupName},
	If[StringLength[$projectPath]>0,
		$fileName=FileNameJoin[{$projectPath,$projectName<>".nb"}];
		dateStamp=DateString[{"YearShort","Month","Day","_","Hour","Minute","Second"}];
		$backupName=FileNameJoin[{$backupPath,$projectName<>"_"<>dateStamp<>".nb"}];
		CopyFile[$fileName,$backupName];
	];
]

ListResources[match_]:=Module[{},
	SetDirectory[$resourcePath];
	Map[{FileNameSplit[#][[-1]],#}&,FileNames[match,"",Infinity]]
];

SaveGraphics[name_,graphic_,opts___]:=Module[{addOptions,standardOptions,graphicsResolution,graphicsList,namesList,graphicsType},
	standardOptions = {"eps"->{},"pdf"->{},"png"->{ImageResolution->150},"tiff"->{"ImageEncoding"->"LZW",ImageResolution->300}};
	graphicsType=GraphicsType/.{opts}/.{GraphicsType->$DefaultGraphicsType};
	graphicsList=If[Head[graphic]=!=List,{graphic},graphic];
	graphicsResolution = Resolution/.({opts}/.{Screen->72, Printer->300, Medium->150, Low->75, High->600})/.{Resolution->$DefaultImageResolution};
	(*Message[SaveGraphics::res, graphicsResolution];*)
	namesList=If[Head[name]=!=List,{name},name];
	addOptions=Options/.{opts}/.Options->{};
	If[Length[namesList]<Length[graphicsList],
		If[Mod[Length[graphicsList],Length[namesList]]==0,
			namesList=MapIndexed[namesList[[1+Mod[#2[[1]],Length[namesList]]]]<>"_"<>StringTake["0000"<>ToString[Floor[(1+#2[[1]])/Length[namesList]]],-4]&,graphicsList];
		]
	];
	If[Length[namesList]>=Length[graphicsList],
		namesList=namesList[[Range[Length[graphicsList]]]];
		Map[
			Table[
				Export[FileNameJoin[{$imagePath,#[[1]]<>"."<>fileType}],#[[2]],fileType,ImageResolution->graphicsResolution,Sequence[fileType/.addOptions/.standardOptions]]
			,{fileType,graphicsType}]&,
			Transpose[{namesList,graphicsList}]
		]
	]
];

SaveAllGraphicsByName[prefix_,saveprefix_]:=Module[{listOfObjects},
	NotebookSave[];
	listOfObjects=Names[prefix<>"*"];
	objects=Map[{saveprefix<>StringTake[#,{StringLength[prefix]+1,-1}],ToExpression[#]}&,listOfObjects];
	Map[SaveGraphics[#[[1]],#[[2]]]&,objects];
]

SaveGraphicsByName[name_,saveprefix_]:=Module[{listOfObjects},
	NotebookSave[];
	listOfObjects={name};
	objects=Map[{saveprefix<>StringTake[#,{StringLength[prefix]+1,-1}],ToExpression[#]}&,listOfObjects];
	Map[SaveGraphics[#[[1]],#[[2]]]&,objects];
]

SetSGContext[]:=Module[{old,new},
	old=$SaveGraphicsContext;
	new=InputString["New Context",old];
	$SaveGraphicsContext=new;
]
SetSGPrefix[]:=Module[{},
	old=$SaveGraphicsPrefix;
	new=InputString["New Prefix",old];
	$SaveGraphicsPrefix=new;
]

MemoryForm[number_]:=Module[{},
	Which[
		number<10^3,ToString[number],
		number<10^6,ToString[Floor[number/1000]]<>"k",
		number<10^9,ToString[Floor[number/10^6]]<>"M"
	]
];



CreateSaveGraphicsPalette[]:=Module[{listOfObjects,content,objects,width=320},
	listOfObjects = Names[$SaveGraphicsContext<>"*"];
	objects = Map[{$SaveGraphicsPrefix<>StringTake[#,{StringLength[$SaveGraphicsContext]+1,-1}],#,StringTake[#,{StringLength[$SaveGraphicsContext]+1,-1}]}&,listOfObjects];
	content = Grid[Join[
			{{Row[{Style["Resolution ",FontFamily->{"Palatino"},FontSize->8],SetterBar[Dynamic[$DefaultImageResolution],Table[res->ToString[res],{res,{75,150,300,600}}],Appearance->"Palette",Alignment->Center,ImageSize->25,BaseStyle->{"Palatino",8}]}]}},
			{{Row[{Style["ImageType ",FontFamily->{"Palatino"},FontSize->8],TogglerBar[Dynamic[$DefaultGraphicsType],Table[res->ToString[res],{res,{"png","pdf","eps","tiff"}}],Appearance->"Palette",Alignment->Center,ImageSize->25,BaseStyle->{"Palatino",8}]}]}},
			{{Row[{Style["FilePrefix ",FontFamily->{"Palatino"},FontSize->8],Button[Dynamic[$SaveGraphicsPrefix],Null,Alignment->Left,Appearance->"Palette",ImageSize->104,BaseStyle->{"Palatino",8}]}]}},
			{{Row[{Style["Context ",FontFamily->{"Palatino"},FontSize->8],Button[Dynamic[$SaveGraphicsContext],Null,Alignment->Left,Appearance->"Palette",ImageSize->104,BaseStyle->{"Palatino",8}]}]}},
			{{Row[{Style["Folder ",FontFamily->{"Palatino"},FontSize->8],
				FileNameSetter[Dynamic[$imagePath],"Directory",ImageSize->{25,18},Appearance->"",BaseStyle->{"Palatino",8}]}]}},
			{{Button[Dynamic[FileNameJoin[FileNameSplit[$imagePath][[{-1}]]]],Null,Alignment->Left,Appearance->"Palette",ImageSize->width,BaseStyle->{"Palatino",8}]}},
			{{Row[{Style["Update ",FontFamily->{"Palatino"},FontSize->8],Button["UPDATE",CreateSaveGraphicsPalette[],Appearance->"Palette",ImageSize->104,BaseStyle->{"Palatino",8,Bold}]}]}},
			{{Row[{Style["",FontFamily->{"Palatino"},FontSize->8],Button["FOLDER",OpenImageFolder[],Appearance->"Palette",ImageSize->104,BaseStyle->{"Palatino",8,Bold}]}]}},
			{{Button["GRAPHICS",Null,Appearance->"Palette",ImageSize->{width,15},BaseStyle->{"Palatino",8,Bold},ContentPadding->False]}},
			Map[{Row[{
				Tooltip[Button[#[[3]],SaveGraphics[#[[1]],ToExpression[#[[2]]]],Alignment->Left,Appearance->"Palette",ImageSize->width-40,BaseStyle->{"Palatino",8},Method->"Queued"],"Save as: "<>#1[[1]]],
				Button[MemoryForm[ByteCount[ToExpression[#[[2]]]]],(Remove[Evaluate[#[[2]]]];CreateSaveGraphicsPalette[];),Alignment->Right,Appearance->"Palette",ImageSize->30, BaseStyle->{"Palatino",8},Method->"Queued"]
			}]}&,objects]
		]
		,Spacings->{0,0},Alignment->Right];
	If[Head[$saveGraphicsPaletteWindow]===Symbol||(Options[$saveGraphicsPaletteWindow,Visible]===$Failed),
		$saveGraphicsPaletteWindow=CreatePalette[content,WindowTitle->"Save Graphics",WindowSize->{320,640}];,
		CreatePalette[content,$saveGraphicsPaletteWindow,WindowTitle->"Save Graphics",WindowSize->{320,640}];
	];
]

ExportProject[opts___]:=Module[{foldersToInclude,folderList,zipCommand,exportName,dateString},
	foldersToInclude=Folders/.{opts}/.Folders->{"images","manuscript"};
	dateString=DateString[{"Year","Month","Day","_","Hour","Minute"}];
	$projectPath=FileNameJoin[Drop[FileNameSplit[NotebookFileName[]],-1]];
	$projectName=FileBaseName[NotebookFileName[]];
	exportName=Name/.{opts}/.Name->$projectName<>dateString;
	folderList=Map[" "<>FileNameJoin[{#}]&,foldersToInclude];
	SetDirectory[$projectPath];
	If[!FileExistsQ["exports"],CreateDirectory[FileNameJoin[{$projectPath,"exports"}]]];
	zipCommand="!zip -r exports"<>$PathnameSeparator<>exportName<>StringJoin[folderList];
	ReadList[zipCommand,String];
	FileByteCount[FileNameJoin[{$projectPath,"exports",exportName<>".zip"}]]
];

OpenProjectFolder[] := Module[{},
	Run["open",$projectPath];
]

OpenImageFolder[] := Module[{},
	Run["open",$imagePath];
]

GetFromAddressBook[name_] := Module[{},
	RunThrough["/opt/local/bin/contacts -H -l -f '%we'", name]
]

CreateAuthor[name_] := Module[{data, contact},
	data = ReadList["!/opt/local/bin/contacts -H -l -f '%fn;%ln;%we;%p;%c, %wa' "<>name,"String"][[1]];
	data=Map[StringTrim,StringSplit[data,";"]];
	contact = Thread[{"Firstname","Lastname","Email","Phone","Address"}->data];
	CellPrint[Cell[("Firstname"/.contact)<>" "<>("Lastname"/.contact),"Author", CellFrameMargins -> 4]];
	CellPrint[Cell[("Address")/.contact,"AuthorAffiliation"]];
	CellPrint[Cell["Email"/.contact,"AuthorEmail"]];
]



End[ ];


EndPackage[ ];
