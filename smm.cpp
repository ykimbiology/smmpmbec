/*
www.mhc-pathway.net/smm
Original file by Bjoern Peters.

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any
damages arising from the use of this software.

Permission is granted to anyone to use this software for any
purpose, including commercial applications, and to alter it and
redistribute it freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must
not claim that you wrote the original software. If you use this
software in a product, an acknowledgment in the product documentation
would be appreciated but is not required.

2. Altered source versions must be plainly marked as such, and
must not be misrepresented as being the original software.

3. This notice may not be removed or altered from any source
distribution.
*/

// 20071115: SMM-BP uses a Bayesian prior derived from randomized peptide library data
//   to model binding specificity of peptide:protein complexes.
//   The program trains a scoring matrix given a list of peptides with corresponding measured
//	 binding affinities.

// Change Log:
// 1.01 changed handling of failed conversion from Exception to ignore if > 2 / 3 of repeat runs converge


#include "stdafx.h"
#include "SMMSeqPair.h"
#include "math.h"

string AutoFileName(const string & original, const string &prefix)
{
	unsigned pos1 = original.find_last_of('\\',original.size());
	unsigned pos2 = original.find_last_of('/',original.size());

	unsigned pos = min(pos1,pos2);
	if(pos==original.size())
		return(prefix + "-" + original);

	if(max(pos1,pos2)<original.size())
		pos=max(pos1,pos2);
	string result;
	result = original.substr(0,pos+1);
	result +=prefix;
	result +="-";
	result +=original.substr(pos+1,original.size()-pos);
	return(string(result));
}

void InvalidParametersMessage(bool usage_only)
{
	if(usage_only)
		cout << "Stablized Matrix Method (SMM-BP) release (2.01)"<< endl;
	else
		cout << "Invalid parameters"<< endl;
	cout << "Usage: smmpmbec [parameter] [input filename] [output filename (optional)]" << endl;
}


void CommandLineInterface(int argc, char* argv[])
{
	cout << endl << "Starting...";
	cout.flush();
	InitParamCV init_param;

	string command=argv[1];
	if(command[0]=='x')
	{
		CNumMat m(20,20);
		m(0,0) = 2.5114372050328058; m(0,1) = -1.172556289749164; m(0,2) = 0.35785719121210213; m(0,3) = -1.0066293422173056; m(0,4) = 0.1727580403054758; m(0,5) = -1.241509446028581; m(0,6) = 0.5054672567192254; m(0,7) = -0.18054599173532332; m(0,8) = 0.5711665055648867; m(0,9) = -0.12111372359481966; m(0,10) = -0.7110238388371026; m(0,11) = 0.8064758371789675; m(0,12) = -0.7111157158338757; m(0,13) = 1.7808611135576138; m(0,14) = -0.8198123695164092; m(0,15) = -1.0429013208622835; m(0,16) = -0.23537449128318744; m(0,17) = -0.38790770600200647; m(0,18) = 2.3586364424295496; m(0,19) = -0.43416935634056786;
		m(1,0) = -1.172556289749164; m(1,1) = 5.161255644919768; m(1,2) = -0.08204952336951452; m(1,3) = -1.218956553546217; m(1,4) = 0.15395860592606125; m(1,5) = -1.914453289989922; m(1,6) = -0.6310519767664662; m(1,7) = -0.6112177498757079; m(1,8) = -0.21103257051666308; m(1,9) = 0.24046091397350455; m(1,10) = 0.5973421443318981; m(1,11) = -1.2198706219703828; m(1,12) = -0.13827959171330167; m(1,13) = 1.4366357124940548; m(1,14) = 0.10154302517352773; m(1,15) = 0.6860602244576559; m(1,16) = 0.19857033555773254; m(1,17) = 1.038438457801885; m(1,18) = -2.3584881241088205; m(1,19) = 0.9436912269700698;
		m(2,0) = 0.35785719121210213; m(2,1) = -0.08204952336951452; m(2,2) = 5.607708036861558; m(2,3) = -3.5456342005577595; m(2,4) = -0.2694064425061858; m(2,5) = 0.9672269599268427; m(2,6) = 0.057798620328543826; m(2,7) = -0.754263146863505; m(2,8) = 0.7004303856692735; m(2,9) = 0.8386737374414306; m(2,10) = -0.07651812326903232; m(2,11) = -0.3891589349000444; m(2,12) = -1.0523488847650042; m(2,13) = 0.5339243740202226; m(2,14) = -0.6682160593383785; m(2,15) = -0.608301257334131; m(2,16) = -0.6598246412078413; m(2,17) = 0.11466364255432451; m(2,18) = 0.7614334304170921; m(2,19) = -0.8339951643199894;
		m(3,0) = -1.0066293422173056; m(3,1) = -1.218956553546217; m(3,2) = -3.5456342005577595; m(3,3) = 5.7053981042845745; m(3,4) = -0.27350259301639257; m(3,5) = 0.9704682509264766; m(3,6) = -0.06064406997015084; m(3,7) = 1.0134567365057752; m(3,8) = -0.24061709488957683; m(3,9) = -1.46968532889567; m(3,10) = 0.3241641724595524; m(3,11) = 0.8432159056369462; m(3,12) = 0.9791098950255378; m(3,13) = -2.3933631849894024; m(3,14) = 1.0511993854426187; m(3,15) = 0.14165271661106532; m(3,16) = 0.7918790712607592; m(3,17) = -0.1373423212069669; m(3,18) = -1.0087425510591945; m(3,19) = 0.5345730021953236;
		m(4,0) = 0.1727580403054758; m(4,1) = 0.15395860592606125; m(4,2) = -0.2694064425061858; m(4,3) = -0.27350259301639257; m(4,4) = 1.521259503903901; m(4,5) = -1.4396251052660958; m(4,6) = -0.3589305312861292; m(4,7) = -0.06300546309897619; m(4,8) = 0.22397035957471556; m(4,9) = -0.5451268912191232; m(4,10) = 0.11615839473140144; m(4,11) = -0.3168596677874526; m(4,12) = 0.1250635892015526; m(4,13) = 1.2249259104507784; m(4,14) = -0.10265622593928761; m(4,15) = 0.6129042680243173; m(4,16) = -0.0005154676473121541; m(4,17) = 0.3988078336611852; m(4,18) = -0.29526748901575317; m(4,19) = 0.11508937100331944;
		m(5,0) = -1.241509446028581; m(5,1) = -1.914453289989922; m(5,2) = 0.9672269599268427; m(5,3) = 0.9704682509264766; m(5,4) = -1.4396251052660958; m(5,5) = 7.222538543399374; m(5,6) = 0.8917540412235465; m(5,7) = 0.03803383447178142; m(5,8) = -0.48137580782017164; m(5,9) = 1.0429715795206707; m(5,10) = 0.34685369870887217; m(5,11) = 0.8508970353345453; m(5,12) = 0.48922353595817436; m(5,13) = -4.032872202024063; m(5,14) = 0.8062086787197222; m(5,15) = -1.6482172765854295; m(5,16) = 0.7218888504076505; m(5,17) = -1.4186742383381525; m(5,18) = 0.08448567207231297; m(5,19) = -1.2558233146175555;
		m(6,0) = 0.5054672567192254; m(6,1) = -0.6310519767664662; m(6,2) = 0.057798620328543826; m(6,3) = -0.06064406997015084; m(6,4) = -0.3589305312861292; m(6,5) = 0.8917540412235465; m(6,6) = 2.2191018685342057; m(6,7) = 0.24045945379642578; m(6,8) = -0.22984526798221777; m(6,9) = -0.18621614852832985; m(6,10) = -0.2876073078433687; m(6,11) = 0.0961795806604124; m(6,12) = 0.04044329212891055; m(6,13) = 0.06244600545387001; m(6,14) = -0.4438142928050018; m(6,15) = -1.043217596169303; m(6,16) = -0.034195403176699286; m(6,17) = 0.23267630402178635; m(6,18) = 0.24525135956039273; m(6,19) = -0.3160551878996512;
		m(7,0) = -0.18054599173532332; m(7,1) = -0.6112177498757079; m(7,2) = -0.754263146863505; m(7,3) = 1.0134567365057752; m(7,4) = -0.06300546309897619; m(7,5) = 0.03803383447178142; m(7,6) = 0.24045945379642578; m(7,7) = 2.7145454086456406; m(7,8) = -0.4076685682500293; m(7,9) = -0.7182685540236934; m(7,10) = 0.025835613101616597; m(7,11) = 0.09753756882233178; m(7,12) = 0.4712110255348285; m(7,13) = -0.4412378899067242; m(7,14) = 0.49223272703468335; m(7,15) = 0.38422738906275716; m(7,16) = 0.4882626901283703; m(7,17) = -1.3149640493529013; m(7,18) = -0.8141224910095863; m(7,19) = 0.3394914570122378;
		m(8,0) = 0.5711665055648867; m(8,1) = -0.21103257051666308; m(8,2) = 0.7004303856692735; m(8,3) = -0.24061709488957683; m(8,4) = 0.22397035957471556; m(8,5) = -0.48137580782017164; m(8,6) = -0.22984526798221777; m(8,7) = -0.4076685682500293; m(8,8) = 2.7317564629981965; m(8,9) = -1.0196026374484826; m(8,10) = -0.2352791846036993; m(8,11) = 0.3129298539261035; m(8,12) = -0.29400602850695023; m(8,13) = 1.1994158774360422; m(8,14) = -1.4952887227418383; m(8,15) = -1.8022213207178028; m(8,16) = 1.1434055715515443; m(8,17) = -0.302894329404214; m(8,18) = 1.0344993078947466; m(8,19) = -0.19774279173386533;
		m(9,0) = -0.12111372359481966; m(9,1) = 0.24046091397350455; m(9,2) = 0.8386737374414306; m(9,3) = -1.46968532889567; m(9,4) = -0.5451268912191232; m(9,5) = 1.0429715795206707; m(9,6) = -0.18621614852832985; m(9,7) = -0.7182685540236934; m(9,8) = -1.0196026374484826; m(9,9) = 5.248888841271222; m(9,10) = -1.2240889829298913; m(9,11) = -0.2821056162349559; m(9,12) = 0.2080241412209852; m(9,13) = -1.7082288821139522; m(9,14) = 0.6127925155499905; m(9,15) = 2.0428778993542798; m(9,16) = -0.8461704937742911; m(9,17) = -1.4608978377636075; m(9,18) = 0.6081006784258332; m(9,19) = -0.26128521023110296;
		m(10,0) = -0.7110238388371026; m(10,1) = 0.5973421443318981; m(10,2) = -0.07651812326903232; m(10,3) = 0.3241641724595524; m(10,4) = 0.11615839473140144; m(10,5) = 0.34685369870887217; m(10,6) = -0.2876073078433687; m(10,7) = 0.025835613101616597; m(10,8) = -0.2352791846036993; m(10,9) = -1.2240889829298913; m(10,10) = 3.453863739264552; m(10,11) = -0.33225867596798286; m(10,12) = 0.7485768400455381; m(10,13) = -1.0228166428064767; m(10,14) = 0.4157512001431434; m(10,15) = -0.2500905756986507; m(10,16) = 0.6365865435654587; m(10,17) = 0.13819556730358307; m(10,18) = -1.5469581495857487; m(10,19) = -0.11668643211366123;
		m(11,0) = 0.8064758371789675; m(11,1) = -1.2198706219703828; m(11,2) = -0.3891589349000444; m(11,3) = 0.8432159056369462; m(11,4) = -0.3168596677874526; m(11,5) = 0.8508970353345453; m(11,6) = 0.0961795806604124; m(11,7) = 0.09753756882233178; m(11,8) = 0.3129298539261035; m(11,9) = -0.2821056162349559; m(11,10) = -0.33225867596798286; m(11,11) = 3.2887321478719875; m(11,12) = 0.1231545309508327; m(11,13) = -1.1089641056258626; m(11,14) = 0.18725982858686246; m(11,15) = -2.3841447893298993; m(11,16) = 1.5383387112722566; m(11,17) = -1.032142805424403; m(11,18) = 0.43933060181183237; m(11,19) = -0.5185463848120979;
		m(12,0) = -0.7111157158338757; m(12,1) = -0.13827959171330167; m(12,2) = -1.0523488847650042; m(12,3) = 0.9791098950255378; m(12,4) = 0.1250635892015526; m(12,5) = 0.48922353595817436; m(12,6) = 0.04044329212891055; m(12,7) = 0.4712110255348285; m(12,8) = -0.29400602850695023; m(12,9) = 0.2080241412209852; m(12,10) = 0.7485768400455381; m(12,11) = 0.1231545309508327; m(12,12) = 1.4998168426878997; m(12,13) = -1.3400565528652237; m(12,14) = 0.41671537036912587; m(12,15) = 0.5179002464455043; m(12,16) = 0.289722943303248; m(12,17) = -0.6927673823839322; m(12,18) = -0.9469309923024158; m(12,19) = 0.26654289549856336;
		m(13,0) = 1.7808611135576138; m(13,1) = 1.4366357124940548; m(13,2) = 0.5339243740202226; m(13,3) = -2.3933631849894024; m(13,4) = 1.2249259104507784; m(13,5) = -4.032872202024063; m(13,6) = 0.06244600545387001; m(13,7) = -0.4412378899067242; m(13,8) = 1.1994158774360422; m(13,9) = -1.7082288821139522; m(13,10) = -1.0228166428064767; m(13,11) = -1.1089641056258626; m(13,12) = -1.3400565528652237; m(13,13) = 6.322943109560347; m(13,14) = -1.6023869169826979; m(13,15) = -0.788163348564028; m(13,16) = -0.6516948621924594; m(13,17) = 1.7400078647832735; m(13,18) = 1.1939258473363354; m(13,19) = 0.5946987729783553;
		m(14,0) = -0.8198123695164092; m(14,1) = 0.10154302517352773; m(14,2) = -0.6682160593383785; m(14,3) = 1.0511993854426187; m(14,4) = -0.10265622593928761; m(14,5) = 0.8062086787197222; m(14,6) = -0.4438142928050018; m(14,7) = 0.49223272703468335; m(14,8) = -1.4952887227418383; m(14,9) = 0.6127925155499905; m(14,10) = 0.4157512001431434; m(14,11) = 0.18725982858686246; m(14,12) = 0.41671537036912587; m(14,13) = -1.6023869169826979; m(14,14) = 2.0426031838217598; m(14,15) = 0.5581502197955588; m(14,16) = 0.2210112409448391; m(14,17) = 0.08200986066375282; m(14,18) = -0.9370738583446581; m(14,19) = 0.08177120942268593;
		m(15,0) = -1.0429013208622835; m(15,1) = 0.6860602244576559; m(15,2) = -0.608301257334131; m(15,3) = 0.14165271661106532; m(15,4) = 0.6129042680243173; m(15,5) = -1.6482172765854295; m(15,6) = -1.043217596169303; m(15,7) = 0.38422738906275716; m(15,8) = -1.8022213207178028; m(15,9) = 2.0428778993542798; m(15,10) = -0.2500905756986507; m(15,11) = -2.3841447893298993; m(15,12) = 0.5179002464455043; m(15,13) = -0.788163348564028; m(15,14) = 0.5581502197955588; m(15,15) = 9.79418781811431; m(15,16) = -6.287946738020439; m(15,17) = 1.641670947597436; m(15,18) = -0.7680252494856169; m(15,19) = 1.2435977433047043;
		m(16,0) = -0.23537449128318744; m(16,1) = 0.19857033555773254; m(16,2) = -0.6598246412078413; m(16,3) = 0.7918790712607592; m(16,4) = -0.0005154676473121541; m(16,5) = 0.7218888504076505; m(16,6) = -0.034195403176699286; m(16,7) = 0.4882626901283703; m(16,8) = 1.1434055715515443; m(16,9) = -0.8461704937742911; m(16,10) = 0.6365865435654587; m(16,11) = 1.5383387112722566; m(16,12) = 0.289722943303248; m(16,13) = -0.6516948621924594; m(16,14) = 0.2210112409448391; m(16,15) = -6.287946738020439; m(16,16) = 7.609681840780615; m(16,17) = -2.9071672734130596; m(16,18) = -0.6809677312286102; m(16,19) = -0.3354906968285788;
		m(17,0) = -0.38790770600200647; m(17,1) = 1.038438457801885; m(17,2) = 0.11466364255432451; m(17,3) = -0.1373423212069669; m(17,4) = 0.3988078336611852; m(17,5) = -1.4186742383381525; m(17,6) = 0.23267630402178635; m(17,7) = -1.3149640493529013; m(17,8) = -0.302894329404214; m(17,9) = -1.4608978377636075; m(17,10) = 0.13819556730358307; m(17,11) = -1.032142805424403; m(17,12) = -0.6927673823839322; m(17,13) = 1.7400078647832735; m(17,14) = 0.08200986066375282; m(17,15) = 1.641670947597436; m(17,16) = -2.9071672734130596; m(17,17) = 5.012199479675745; m(17,18) = -0.40970664183537625; m(17,19) = 0.6657946270616542;
		m(18,0) = 2.3586364424295496; m(18,1) = -2.3584881241088205; m(18,2) = 0.7614334304170921; m(18,3) = -1.0087425510591945; m(18,4) = -0.29526748901575317; m(18,5) = 0.08448567207231297; m(18,6) = 0.24525135956039273; m(18,7) = -0.8141224910095863; m(18,8) = 1.0344993078947466; m(18,9) = 0.6081006784258332; m(18,10) = -1.5469581495857487; m(18,11) = 0.43933060181183237; m(18,12) = -0.9469309923024158; m(18,13) = 1.1939258473363354; m(18,14) = -0.9370738583446581; m(18,15) = -0.7680252494856169; m(18,16) = -0.6809677312286102; m(18,17) = -0.40970664183537625; m(18,18) = 5.479998150959398; m(18,19) = -1.4393782129317112;
		m(19,0) = -0.43416935634056786; m(19,1) = 0.9436912269700698; m(19,2) = -0.8339951643199894; m(19,3) = 0.5345730021953236; m(19,4) = 0.11508937100331944; m(19,5) = -1.2558233146175555; m(19,6) = -0.3160551878996512; m(19,7) = 0.3394914570122378; m(19,8) = -0.19774279173386533; m(19,9) = -0.26128521023110296; m(19,10) = -0.11668643211366123; m(19,11) = -0.5185463848120979; m(19,12) = 0.26654289549856336; m(19,13) = 0.5946987729783553; m(19,14) = 0.08177120942268593; m(19,15) = 1.2435977433047043; m(19,16) = -0.3354906968285788; m(19,17) = 0.6657946270616542; m(19,18) = -1.4393782129317112; m(19,19) = 1.9239224463818685;
		init_param.inverse_covar = m;
		command=command.substr(1,command.size()-1);
		cout << command;
	}
	if(command=="-cone" || command=="-cgroup" || command=="-lone" || command=="-lgroup" || command=="-coneEqual")
	{
		bool allEqual = false;	// false = inequalities are kept; true = all signs are set to equal.

		if( command == "-cone")
			init_param.lambda_grouping=ONE_COVAR;
		else if(command == "-cgroup")
			init_param.lambda_grouping=GROUP_COVAR;
		else if(command == "-lone")
			init_param.lambda_grouping=ONE_LAMBDA;
		else if (command == "-lgroup")
			init_param.lambda_grouping=GROUP_LAMBDA;
		else if (command == "-coneEqual") {
			init_param.lambda_grouping=ONE_COVAR;
			allEqual = true;
		}


		const string in_file_name = argv[2];	// the input file containing the binding data for peptides.
		CSeqSet set;
		set.Load(in_file_name, allEqual);
		CSeqMatrix mat;
		mat.SMMTrain(set,init_param);

		// DEBUG:
		// cout << "mat size " << mat.GetMatrixSize() << endl;

		if (mat.GetMatrixSize() > 0) {
			if(argc < 4)
				mat.Save(AutoFileName(in_file_name,"mat"));
			else
				mat.Save(argv[3]);
		}

	}
	else if (command=="-production") {
		// This is the version that a user will run.
		// Here is the order of attempts:
		// 1. -cgroup
		// 2. -cone
		// 3. -cone with all inequalities changed to equal signs.
		// Step 1 will be tried first and remaining steps will be executed in order
		//  if the previous step does not yield a trained scoring matrix and throws an exception.

		const string in_file_name=argv[2];	// the input file containing the binding data for peptides.
		CSeqSet set;

		bool allEqual = false;
		init_param.lambda_grouping=GROUP_COVAR;
		set.Load(in_file_name, allEqual);
		CSeqMatrix mat;


		try {
			mat.SMMTrain(set,init_param);
		} catch (BPNoConvergence e) {

			// If first try with -cgroup option fails, write out a log file.
			string output_name = "log-file.txt";
			ex_ofstream logfile(output_name);
			logfile << "clog: -production" << endl << endl;
			logfile << "1st attempt (-cgroup) failed." << endl << endl;
			logfile << e.m_message.str() << endl << endl;

			try {
				//DEBUG:
				logfile << "Trying 2nd attempt (-cone) ..." << endl;
				logfile << "allEqual "<< allEqual << endl;
				logfile << "matrix size "<< mat.GetMatrixSize() << endl;
				logfile << "init_param: " << init_param.lambda_grouping << endl << endl;

				init_param.lambda_grouping=ONE_COVAR;
				set.Load(in_file_name, allEqual);
				mat.SMMTrain(set,init_param);
			} catch (BPNoConvergence e) {
				//DEBUG:
				logfile << "2nd attempt (-cone) failed." << endl;
				logfile << "Trying 3rd attempt (-cone + conversion of all inequal signs to equal signs) ..." << endl;
				logfile << "allEqual "<< allEqual << endl;
				logfile << "matrix size "<< mat.GetMatrixSize() << endl;
				logfile << "init_param: " << init_param.lambda_grouping << endl;

				init_param.lambda_grouping=ONE_COVAR;
				allEqual = true;	// inequalities are converted into equal signs.
				set.Load(in_file_name, allEqual);
				mat.SMMTrain(set,init_param);
			}
		}

		// If matrix has been trained, then write it out.
		if (mat.GetMatrixSize() > 0) {
			if(argc <4)
				mat.Save(AutoFileName(in_file_name,"mat"));
			else
				mat.Save(argv[3]);
		}
	}
	else
		throw BPException(string("Unknown command line parameter '"+command+"'"));

	cout << endl << "Executed succesfully" << endl;
}



void Background()
{
	ex_chdir("Background");

	bool allEqual = false;
	CSeqSet set;
	set.Load("RawData/murine-tap-9-mers_NO_EXTREME_OUTLIERS_NO_INEQUAL.txt", allEqual);

	InitParamCV init_param;
	init_param.lambda_grouping=ONE_COVAR;

	CSMMSet smmset;
	set.ConvertToSMMSet(smmset, true);

	CSMMSolve smmsolve;
	smmsolve.InitSolver(smmset, init_param);
	CNumVec param(1);
	param[0]=0;
	smmsolve.SolveX(smmset, param);
	clog << endl << "Solution for entire set with Param = 0";
	clog << endl << smmsolve.GetLambda();
	clog << endl << smmsolve.GetX();

	param[0]=-1;
	smmsolve.SolveX(smmset, param);
	clog << endl << "Solution for entire set with Param = -1";
	clog << endl << smmsolve.GetLambda();
	clog << endl << smmsolve.GetX();

	param[0]=1;
	smmsolve.SolveX(smmset, param);
	clog << endl << "Solution for entire set with Param = 1";
	clog << endl << smmsolve.GetLambda();
	clog << endl << smmsolve.GetX();
	// Confirmed with Mathematica


	clog << endl << endl << "Next, compare distance in cross validation. For this, shuffling in CVSet has to be disabled. Don't forget to re-enable";
	init_param.lambda_grouping=ONE_COVAR;
	init_param.cv_num=5;
	init_param.train_repeats=1;
	CSeqMatrix mat;
	mat.SMMTrain(set,init_param);
	mat.Save("mat-ONE_COVAR.txt");
	// Confirmed with Mathematica that this is minimum

	clog << endl << endl << "Compare the Gradient at All Lambda = 0.264305";
	init_param.lambda_grouping=GROUP_COVAR;
	init_param.precision=0.0000001;
	mat.SMMTrain(set,init_param);
	mat.Save("mat-GROUP_COVAR.txt");
	// confirmed with Mathematica that this is minimum
}


int main(int argc, char* argv[])
{
	try
	{
		// Redirect clog, cerr to strstreams, which can later be inserted into the output file
		stringstream output_clog, output_cerr;
		streambuf* clog_old = clog.rdbuf(output_clog.rdbuf());
		streambuf* cerr_old = cerr.rdbuf(output_cerr.rdbuf());
		// Default output detail
		clog_detail.push_back(MEDIUM);
		string output_name = "log-file.txt";

		CLetter::Init("ACDEFGHIKLMNPQRSTVWY");	// Assume peptide sequences as input

		try
		{
			if(argc < 2) // Catch trivial missing parameters --> Display usage information
			{
				InvalidParametersMessage(true);
				return(0);
			}
			string command = argv[1];
			if(command[0]!='-' && command[0]!='x')
				InvalidParametersMessage(false);
			else if(command=="-background")
			{
				output_name = "CPP-background.txt";
				Background();
				ex_ofstream logfile(output_name);
				logfile << "SelfTest Commpleted without Exception" << endl;
				logfile << "CLOG:" << endl;
				logfile << output_clog.str();
				logfile << endl << "CERR" << endl;
				logfile << output_cerr.str();
			}
			else if(argc<3)
				InvalidParametersMessage(false);
			else
				CommandLineInterface(argc,argv);
			clog.rdbuf(clog_old);
			cerr.rdbuf(cerr_old);
			return(0);
		}
		catch(BPFileException e)
		{
		    cerr << endl<< endl<< "File related exception, program terminated";
		    cerr << endl<< endl<< "Message: " <<  e.m_message.str();
		    cerr << endl<< endl<< "File: " << e.m_filename;
		    cerr << endl<< endl<< "Current Directory: " << ex_getcwd();
		}
		catch(BPException e)
		{
		    e.m_message << endl << '\0';
		    cerr << endl<< endl<< "Exception, program terminated";
		    cerr << endl<< endl<< "Message:" <<  e.m_message.str();
		}
		catch(...)
		{
		    cerr << endl<< "Some other exception.";
		}
		cout << endl;
		cout << "An error occured:" << endl;
		cout << output_cerr.str() << endl;

		ex_ofstream logfile(output_name);
		logfile << "clog: inside the outer try block" << endl;
		logfile << output_clog.str();
		logfile << endl << "cerr" << endl;
		logfile << output_cerr.str();

		clog.rdbuf(clog_old);
		cerr.rdbuf(cerr_old);
		return(-1);
	}
	catch(...)
	{
		cerr << endl<< "Exception. I am really not supposed to occur. cerr";
		cout << endl<< "Exception. I am really not supposed to occur. cout";
		return (-2);
	}
	cerr << endl << "Now this should never happen.";
	cout << endl << "Now this should never happen.";
	return(-3);
}



