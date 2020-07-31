# Dilatometry-analysis
The Master Fitter code can be used to analyzed dilatometry results on hypoeutectoid steels. This code uses a technique called Unit Cell Dilation (UCD). For detailed explanation of UCD please refer to the paper <a href="https://www.sciencedirect.com/science/article/abs/pii/S1044580317324051">"Unit cell dilation technique for analyzing dilatometry data in microalloyed steels"</a>. For even more in depth explanation of the method please refer the my PhD thesis titled <a href="https://era.library.ualberta.ca/items/079669b6-d108-4c89-86ed-4e067600c57c">"A fundamental method to quantify phase change in microalloyed steels"</a>. 

<h1> What can Master Fitter do? </h1>
Master Fitter can analyze a series of continuous cooling dilatometry data and calculate phase fractions by curve fitting. To do so, Master Fitter calculates some material properties that are used to calculate amount of sample dilation as a function of amount austenite that is transformed. The product of transformation is determined based on thermodynamic model or empirical equation. An example of thermodynamic model is triggering Fe3C formation when austenite’s carbon content goes above solubility limit of austenite with respect to carbon. An example of empirical equation is triggering bainitic ferrite formation once bainite start temperature is reached.


<h1> Dilatometry basics </h1>
Dilatometry is a material characterization technique in which change in dimension of a sample is recorded as a function of temperature and time. It is one the most accurate and sensitive methods to study phase transformation in materials especially when there is a relatively large difference in the density of the phases. One such material is steel in which a high density austenite phase transforms to multitude of lower density microconstituents and phases such polygonal ferrite, pearlite, bainitic ferrite, bainite, martensite and retained austenite. In simple phase transformation interpreting dilatometry data is relatively easy using so called "lever rule". However, for a complex system such as microallyed steel which often produces more than two kinds of product phase "lever rule" is out right wrong. The reason for this is explained in <a href="https://www.sciencedirect.com/science/article/abs/pii/S1044580317324051">"Unit cell dilation technique for analyzing dilatometry data in microalloyed steels"</a> and in page 82 of <a href="https://era.library.ualberta.ca/items/079669b6-d108-4c89-86ed-4e067600c57c">"A fundamental method to quantify phase change in microalloyed steels"</a>.


<h1> Using Master Fitter </h1>
<body lang="en-CA" link="#000080" vlink="#800000" dir="ltr"><p style="margin-bottom: 0cm; line-height: 100%">
The Master Fitter program is the optimization code that performs the curve fitting using the Unit Cell Dilation algorithm. It can be used in two modes with optimization being on or off. In the optimizationoff mode, Master Fitter just performs curve fitting using the predefined materials properties. In the optimization off mode the output of the Master Fitter is the phase fraction of the phases that form during cooling. In the optimization on mode, Master Fitter searches for the best combination of materials properties that results in the best curve fitting at the same time for a number of dilatometry data files. In this mode the results of Master Fitter are optimized material properties as well as phase fractions of the phases during transformation for all data files.</p>
<p style="margin-bottom: 0cm; line-height: 100%">Master Fitter uses several input files that need to be provided by the user. These files are:</p>
<ul>
	<li><p style="margin-bottom: 0.35cm; line-height: 100%"><b>Dilation
	data files:</b> These files contain dilation data in *.csv format.
	It is highly recommended that these files be names in the specific
	format that tells the program about the nature of the experiment and
	provides some key information. So far the only type of experiment
	that can be analyzed is continuous cooling. Continuous cooling data
	files names should begin with CC and followed by the cooling rate.
	For example <i>CC5.csv</i> is continuous cooling experiment at the
	cooling rate of 5<font face="Times New Roman, serif">°C</font>/s.
	Cooling rates less than 1 are written without the decimal point. For
	example <i>CC02.csv</i> means cooling rate at 0.2<font face="Times New Roman, serif">°C</font>/s
	and <i>CC0004.csv</i> means 0.004<font face="Times New Roman, serif">°C</font>/s.</p>
	<li><p style="margin-bottom: 0.35cm; line-height: 100%"><b>Chemistry</b>:
	Chemical composition of steel is written in the text file named
	<i>chemistry.txt</i>. In this file each line contains the
	concentration of one element in the format <i>element,wt%,</i>. The
	commas are important and needed. In future versions this file format
	can be replaced with .csv or .xlsx for simplicity and consistency.
	The elements that can be defined are C, Mn, Si, Ni, Cr, Mo, Al, Nb,
	Ti, Al and Cu.</p>
	<li><p style="margin-bottom: 0.35cm; line-height: 100%"><b>Carbon
	solubility in austenite</b>: This file contains the solubility data
	used for determining whether carbon concentration in austenite
	permits precipitation of cementite. The file name should be
	<i>C_solubility_limit_austenite.xlsx</i>. Data is stored in the form
	of a Microsoft Excel file with the first column being carbon mole
	fraction and the second column being the temperature in <font face="Times New Roman, serif">°C</font>.</p>
	<li><p style="margin-bottom: 0.35cm; line-height: 100%"><b>Master
	Fitter setup file: </b>This file contains the specific information
	about each data file that Master Fitter needs to understand the
	data. This file is named <i>master_setup.txt</i>. In this file each
	line contains the setting for one data file written in the following
	format:</p>
</ul>
<p style="margin-left: 2.03cm; margin-bottom: 0.35cm; line-height: 100%">
<b>File name, full austenite microstructure start row number, full austenite microstructure finish row number, full product microstructure start row number, full product microstructure finish row number, sample length in meter, cementite instructions</b></p>
<p style="margin-left: 2.03cm; margin-bottom: 0.35cm; line-height: 100%">
<br/>

</p>
<p style="margin-left: 2.03cm; margin-bottom: 0.35cm; line-height: 100%">
The file name must match an existing dilation data file stored in the
same folder as the Master Fitter. Row numbers are different than the
row numbers found in the original dilation data files. Row numbers
should be determined in the plots that <i>data_plotter_Vx.py</i>
generates. In these plots row numbers of the points are shown by each
point. This is because of the conditioning that
<i>data_set_conditioner_Vx.py</i> does to the original data. Sample
length is the length over which dilation is measured. If dilation is
measured over the diameter of the specimen, the diameter is provided
and if dilation is measured over length then specimen length is
provided. Cementite instruction tells the program to allow, disallow
(e.g., high Si) precipitation of cementite. It is also possible to
determine a temperature below which cementite is not allowed to form.
The syntax for these three conditions are:</p>
<ul>
	<ul>
		<li><p style="margin-bottom: 0.35cm; line-height: 100%">Cementite_yes</p>
		<li><p style="margin-bottom: 0.35cm; line-height: 100%">Cementite_no</p>
		<li><p style="margin-bottom: 0.35cm; line-height: 100%">User_XXX
		(XXX is the lowest temperature for cementite precipitation in <font face="Times New Roman, serif">°C</font>)</p>
	</ul>
</ul>
<p style="margin-left: 2.03cm; margin-bottom: 0.35cm; line-height: 100%">
<i>Cementite_yes</i> should be generally be used. The other two conditions must be used only if the user is certain they apply.</p>
