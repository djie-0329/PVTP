#include "stdafx.h"
#include "geodesic_mesh.h"
#include "geodesic_algorithm_exact.h" 
#include "geodesic_algorithm_parallel_fwp_exact.h"
#include "geodesic_algorithm_approximate.h"
#include "geodesic_algorithm_parallel_fwp_approximate.h"
#include "geodesic_result_process.h"

using namespace std;

unsigned Kmin = 6, Kmax = 1e9;
double step = 1.0;
const double eta = 100; //?
double tau_value = 1.0;

//char file_name[255] = { '\0' };

bool NeedTexture = false; //-t
bool NeedRelativeErrorTexture = false; //-r
bool NeedJudge = false; //-judge
bool NeedOutputResult = false; //-o
bool NeedComputeRelativeError = false; //-c

unsigned alg_id = 10;
int source_id = -1;
unsigned N_source = 0;
double lambda = 128.0;
unsigned num_procs = 1;
unsigned K_concurrent = 100;


std::string inFile, outResult;
std::string sourceFile;
geodesic::Mesh mesh;
std::vector<int> sources_list;

void ComputeTau()
{

	double avg_tau = 0, max_tau = 0, var_tau = 0;
	double lenA, lenB, lenC;
	double r, R; //r: in_circle; R: circumradius
	double S, H;//S: area; H: max length
	double p, tau, dau;
	std::vector<double> Taus;
	Taus.clear();

	double avg_dau = 0, max_dau = 0, var_dau = 0;
	std::vector<double> Daus;
	Daus.clear();

	double min_len = 1e12;
	for (unsigned i = 0; i < mesh.edges().size(); i++)
	{
		double len = mesh.edges()[i].length();
		min_len = min(min_len, len);
	}

	for (unsigned i = 0; i < mesh.faces().size(); i++)
	{
		lenA = mesh.faces()[i].adjacent_edges()[0]->length();
		lenB = mesh.faces()[i].adjacent_edges()[1]->length();
		lenC = mesh.faces()[i].adjacent_edges()[2]->length();
		H = 1e-12;
		H = max(lenA, H);
		H = max(lenB, H);
		H = max(lenC, H);
		p = (lenA + lenB + lenC) / 2.0; //semi-perimeter
		S = sqrt(p*(p - lenA)*(p - lenB)*(p - lenC));
		r = sqrt((p - lenA)*(p - lenB)*(p - lenC) / p);
		//tau = 6.0 / sqrt(3.0) * r / max(max(lenA, lenB), lenC);
		//tau = 1.0 / tau;
		tau = (p*H) / (2 * sqrt(3.0)*S);
		dau = tau * max(max(lenA, lenB), lenC) / min_len;
		avg_tau += tau;
		avg_dau += dau;
		max_tau = max(max_tau, tau);
		max_dau = max(max_dau, dau);
		Taus.push_back(tau);
		Daus.push_back(dau);
	}
	avg_tau /= mesh.faces().size();
	avg_dau /= mesh.faces().size();
	for (unsigned i = 0; i < mesh.faces().size(); i++)
	{
		var_tau += (Taus[i] - avg_tau)*(Taus[i] - avg_tau);
		var_dau += (Daus[i] - avg_dau)*(Daus[i] - avg_dau);
	}
	var_tau /= mesh.faces().size();
	var_dau /= mesh.faces().size();
	tau_value = double(int(avg_tau));

	printf("avg_tau: %.6lf, max_tau: %.6lf, var_tau: %.6lf\n", avg_tau, max_tau, var_tau);
}

bool InputParameter(int argc, char **argv)
{
	if (argc == 1)
	{
		cout << endl << "OPTIONS: " << endl;
		cout << endl;

		cout << "-alg [algorithm]: input the index of algorithm. " << endl << endl;
		cout << "VTP: -alg 0 -m -s(-stxt) (-o) (-t);" << endl;
		cout << "Parallel_FWP_VTP: -alg 1 -m -s(-stxt) -np -k (-o) (-t) (-judge);" << endl;
		cout << "Appr_VTP: -alg 2 -m -s(-stxt) -l (-o) (-t) (-c) (-r);" << endl;
		cout << "Parallel_FWP_Appr_VTP: -alg 3 -m -s(-stxt) -l -np -k (-o) (-t) (-c) (-r)." << endl << endl;


		cout << "-m [meshFile]: input model file." << endl << endl;
		cout << "-s [src]: index of source." << endl << endl;
		cout << "-stxt [src txt]: the text of sources index." << endl << endl;
		cout << "-l [lambda]: parameter for Approximate VTP algorithm, repetitive gap" << endl << endl;
		cout << "-np [num_procs]: number of process" << endl << endl;
		cout << "-o [output]: bool: you need to output geodesic distance result" << endl << endl;
		cout << "-t [texture]: bool: you need to output texture file for geodesic distance" << endl << endl;
		cout << "-c [to compute relatvie error]: bool: you need to compute the relative error for each vertex and average" << endl << endl;
		cout << "-r [relative error texture]: bool: you need to output texture file for relative error of approximate geodesic result" << endl << endl;
		cout << "-judge [to judge is exact or not]: bool: you need to judge the result is exact or not" << endl << endl;

		return false;
	}

	for (int i = 1; i < argc;)
	{
		printf("%s\n", argv[i]);
		if (strcmp(argv[i], "-alg") == 0)
		{
			alg_id = atoi(argv[i + 1]);
			printf("%s\n", argv[i+1]);
			i += 2;
		}
		else if (strcmp(argv[i], "-m") == 0)
		{
			inFile = argv[i + 1];
			//strcpy(file_name, argv[i + 1]);
			printf("%s\n", argv[i + 1]);
			i += 2;
		}
		else if (strcmp(argv[i], "-s") == 0)
		{
			source_id = atoi(argv[i + 1]);
			printf("%s\n", argv[i + 1]);
			i += 2;
		}
		else if (strcmp(argv[i], "-stxt") == 0)
		{
			std::string source_str = "_" + std::to_string(N_source) + "sample.txt";
			sourceFile = inFile.substr(0, inFile.length() - 4);
			N_source = atoi(argv[i + 1]);
			printf("%s\n", argv[i + 1]);
			sourceFile.append("_" + std::to_string(N_source) + "sample.txt");
			i += 2;
		}
		else if (strcmp(argv[i], "-l") == 0)
		{
			lambda = atof(argv[i + 1]);
			printf("%s\n", argv[i + 1]);
			i += 2;
		}
		else if (strcmp(argv[i], "-np") == 0)
		{
			num_procs = atoi(argv[i + 1]);
			printf("%s\n", argv[i + 1]);
			K_concurrent = num_procs * 100;
			i += 2;
		}
		else if (strcmp(argv[i], "-o") == 0)
		{
			NeedOutputResult = true;
			i++;
		}
		else if (strcmp(argv[i], "-t") == 0)
		{
			NeedTexture = true;
			i++;
		}
		else if (strcmp(argv[i], "-c") == 0)
		{
			NeedComputeRelativeError = true;
			i++;
		}
		else if (strcmp(argv[i], "-r") == 0)
		{
			NeedRelativeErrorTexture = true;
			i++;
		}
		else if (strcmp(argv[i], "-judge") == 0)
		{
			NeedJudge = true;
			i++;
		}
		else ++i;
	}

	return true;
}

//-alg 0 -m -s (-o) (-t);
void RunVTP()
{
	geodesic::GeodesicAlgorithmExact algorithmVTP(&mesh);

	for (unsigned i = 0; i < sources_list.size(); i++)
	{
		//propagate
		algorithmVTP.propagate(sources_list[i]);

		//-o: output geodesic distances
		if (NeedOutputResult)
		{
			std::string source_str = "_source" + std::to_string(sources_list[i]);
			outResult = inFile.substr(0, inFile.length() - 4);
			outResult.append(source_str);
			outResult.append("_VTP_result.txt");
			ofstream fout(outResult.c_str());
			for (unsigned j = 0; j < mesh.vertices().size(); j++)
			{
				double distance = mesh.vertices()[j].geodesic_distance();
				fout << setprecision(20) << distance << endl;
			}
			fout << endl;
			fout.close();

			//-t: need texture
			if (NeedTexture)
			{
				GeodesicResultProcess ResultProcess;
				ResultProcess.GeodesicResultTexture(inFile, outResult);
			}
		}

		//output data to excel file
		FILE *file = fopen("VTP.csv", "a");
		fprintf(file, "%s, %d, ", inFile.c_str(), sources_list[i]);
		fclose(file);

		algorithmVTP.print_statistics();
	}
	return;
}

//-alg 1 -m -s -np (-o) (-t) (-judge);
void RunParallel_VTP()
{

	tbb::task_scheduler_init init(num_procs); 

	geodesic::GeodesicAlgorithmParallelFWPExact algorithmParallel_FWP_VTP(&mesh, num_procs, K_concurrent);

	// Apply for the 'bucket' structure of FWP 
	algorithmParallel_FWP_VTP.Kmin = Kmin;
	algorithmParallel_FWP_VTP.Kmax = Kmax;
	algorithmParallel_FWP_VTP.step = step;
	algorithmParallel_FWP_VTP.binWidth = mesh.avg_edge()/tau_value; // sqrt(mesh.vertices().size());

	for (unsigned i = 0; i < sources_list.size(); i++)
	{
		//propagate
		algorithmParallel_FWP_VTP.propagate(sources_list[i]);

		//-o: output geodesic distance
		if (NeedOutputResult)
		{
			std::string source_str = "_source" + std::to_string(sources_list[i]);
			outResult = inFile.substr(0, inFile.length() - 4);
			outResult.append(source_str);
			outResult.append("_Parallel_FWP_VTP_result.txt");
			ofstream fout(outResult.c_str());
			for (unsigned j = 0; j < mesh.vertices().size(); j++)
			{
				double distance = mesh.vertices()[j].geodesic_distance();
				fout << setprecision(20) << distance << endl;
			}
			fout << endl;
			fout.close();

			//-t: need texture
			if (NeedTexture)
			{
				GeodesicResultProcess ResultProcess;
				ResultProcess.GeodesicResultTexture(inFile, outResult);
			}

			//-judge: need judge if is exact or not
			if (NeedJudge)
			{
				std::string VTPresult = inFile.substr(0, inFile.length() - 4);
				VTPresult.append(source_str);
				VTPresult.append("_VTP_result.txt");
				GeodesicResultProcess ResultProcess;
				bool isExact = ResultProcess.GeodesicDistanceResultExact(VTPresult, outResult);
				if (!isExact)
				{
					printf("Wrong answer!\n");
					system("pause");
				}
			}
		}

		//output data to excel file
		FILE *file = fopen("Parallel_FWP_VTP.csv", "a");
		fprintf(file, "%s, %d, %d, ", inFile.c_str(), sources_list[i], num_procs);
		fclose(file);

		algorithmParallel_FWP_VTP.print_statistics();
	}
	return;
}

//-alg 2 -m -s -l (-o) (-t) (-c) (-r);
void RunAppr_VTP()
{
	geodesic::GeodesicAlgorithmApproximate algorithmAppr_VTP(&mesh, lambda);
	//printf("lambda: %.3lf\n", lambda);
	for (unsigned i = 0; i < sources_list.size(); i++)
	{
		//propagate
		algorithmAppr_VTP.propagate(sources_list[i]);

		double avg_error_percent = 0;
		double mx_relative_error_percent = 0;
		//-o: output geodesic distances
		if (NeedOutputResult)
		{
			std::string source_str = "_source" + std::to_string(sources_list[i]);
			std::string lambda_str = "_lambda" + std::to_string(lambda);
			outResult = inFile.substr(0, inFile.length() - 4);
			outResult.append(source_str);
			outResult.append(lambda_str);
			outResult.append("_Appr_VTP_result.txt");
			ofstream fout(outResult.c_str());
			for (unsigned j = 0; j < mesh.vertices().size(); j++)
			{
				double distance = mesh.vertices()[j].geodesic_distance();
				fout << setprecision(20) << distance << endl;
			}
			fout << endl;
			fout.close();

			GeodesicResultProcess ResultProcess;

			//-t: need texture
			if (NeedTexture)
			{
				ResultProcess.GeodesicResultTexture(inFile, outResult);
			}

			//-c: compute relative error value
			if (NeedComputeRelativeError)
			{
				std::string VTPresult = inFile.substr(0, inFile.length() - 4);
				VTPresult.append(source_str);
				VTPresult.append("_VTP_result.txt");
				std::string errorFile = outResult.substr(0, outResult.length() - 4);
				errorFile.append("_relativeErrorValue.txt");

				avg_error_percent = ResultProcess.GeodesicDistanceApproximateError(VTPresult, outResult, errorFile);

				//-r: need do relative error value texture
				if (NeedRelativeErrorTexture)
				{
					mx_relative_error_percent = ResultProcess.GeodesicResultRelativeErrorValueTexture(inFile, errorFile);
				}
			}

		}
		else// straightly compute average relative error
		{
			std::string VTPresult = inFile.substr(0, inFile.length() - 4);
			std::string source_str = "_source" + std::to_string(sources_list[i]);
			VTPresult.append(source_str);
			VTPresult.append("_VTP_result.txt");
			std::ifstream fin(VTPresult);
			if (!fin)
			{
				printf("Error: unable to open exact result file: %s\n", VTPresult.c_str());
				//system("pause");
			}
			double res;
			double eps = 1e-12;
			double error;
			unsigned cn_v = 0;
			while (fin >> res)
			{
				if (res < eps)
				{
					cn_v++;
					continue;
				}
				error = (mesh.vertices()[cn_v].geodesic_distance() - res) / res;
				error *= 100;
				if (error > mx_relative_error_percent)
					mx_relative_error_percent = error;
				avg_error_percent += error;
				cn_v++;
			}
			fin.close();
			avg_error_percent /= double(cn_v);
			printf("Relative Error: %.6lf%%\n", avg_error_percent);
		}

		//output data to excel file
		FILE *file = fopen("Appr_VTP.csv", "a");
		fprintf(file, "%s, %d, %d, %.6lf, %.6lf, ", inFile.c_str(), sources_list[i], lambda, avg_error_percent, mx_relative_error_percent);
		fclose(file);

		algorithmAppr_VTP.print_statistics();
	}

	return;
}

//-alg 3 -m -s -l -np (-o) (-t) (-c) (-r)
void RunParallel_Appr_VTP()
{
	tbb::task_scheduler_init init(num_procs);
	geodesic::GeodesicAlgorithmParallelFWPApproximate algorithmParallel_FWP_Appr_VTP(&mesh, lambda, num_procs, K_concurrent);

	// Apply for the 'bucket' structure of FWP 
	algorithmParallel_FWP_Appr_VTP.Kmin = Kmin;
	algorithmParallel_FWP_Appr_VTP.Kmax = Kmax;
	algorithmParallel_FWP_Appr_VTP.step = step;
	algorithmParallel_FWP_Appr_VTP.binWidth = mesh.avg_edge()/tau_value;

	for (unsigned i = 0; i < sources_list.size(); i++)
	{
		//propagate
		algorithmParallel_FWP_Appr_VTP.propagate(sources_list[i]);

		double avg_error_percent = 0;
		double mx_relative_error_percent = 0;

		//-o: output geodesic distance
		if (NeedOutputResult)
		{
			std::string source_str = "_source" + std::to_string(sources_list[i]);
			std::string lambda_str = "_lambda" + std::to_string(lambda);
			outResult = inFile.substr(0, inFile.length() - 4);
			outResult.append(source_str);
			outResult.append(lambda_str);
			outResult.append("_Parallel_FWP_Appr_VTP_result.txt");
			ofstream fout(outResult.c_str());
			for (unsigned j = 0; j < mesh.vertices().size(); j++)
			{
				double distance = mesh.vertices()[j].geodesic_distance();
				fout << setprecision(20) << distance << endl;
			}
			fout << endl;
			fout.close();

			GeodesicResultProcess ResultProcess;
			//-t: need texture
			if (NeedTexture)
			{
				ResultProcess.GeodesicResultTexture(inFile, outResult);
			}

			//-c: compute relative error value
			if (NeedComputeRelativeError)
			{
				std::string VTPresult = inFile.substr(0, inFile.length() - 4);
				VTPresult.append(source_str);
				VTPresult.append("_VTP_result.txt");
				std::string errorFile = outResult.substr(0, outResult.length() - 4);
				errorFile.append("_relativeErrorValue.txt");

				avg_error_percent = ResultProcess.GeodesicDistanceApproximateError(VTPresult, outResult, errorFile);

				//-r: need do relative error value texture
				if (NeedRelativeErrorTexture)
				{
					mx_relative_error_percent = ResultProcess.GeodesicResultRelativeErrorValueTexture(inFile, errorFile);
				}
			}
		}
		else // straightly compute average relative error
		{
			std::string VTPresult = inFile.substr(0, inFile.length() - 4);
			std::string source_str = "_source" + std::to_string(sources_list[i]);
			VTPresult.append(source_str);
			VTPresult.append("_VTP_result.txt");
			std::ifstream fin(VTPresult);
			if (!fin)
			{
				printf("Error: unable to open exact result file: %s\n", VTPresult.c_str());
				//system("pause");
			}
			double res;
			double eps = 1e-12;
			double error;
			unsigned cn_v = 0;
			while (fin >> res)
			{
				if (res < eps)
				{
					cn_v++;
					continue;
				}
				error = (mesh.vertices()[cn_v].geodesic_distance() - res) / res;
				error *= 100;
				if (error > mx_relative_error_percent)
					mx_relative_error_percent = error;
				avg_error_percent += error;
				cn_v++;
			}
			fin.close();
			avg_error_percent /= double(cn_v);
			printf("Relative Error: %.6lf%%\n", avg_error_percent);
		}

		//output data to excel file
		FILE *file = fopen("Parallel_FWP_Appr_VTP.csv", "a");
		fprintf(file, "%s, %d, %d, %d, %.6lf, %.6lf, ", inFile.c_str(), sources_list[i], num_procs, lambda, avg_error_percent, mx_relative_error_percent);
		fclose(file);

		algorithmParallel_FWP_Appr_VTP.print_statistics();
	}
	return;
}


void RunAlgorithmByID(int arg_id)
{
	//multi sources
	if (source_id < 0)
	{
		std::ifstream fin(sourceFile);
		if (!fin)
		{
			sources_list.clear();
			if (N_source > mesh.vertices().size())
			{
				for (unsigned i = 0; i < mesh.vertices().size(); i++)
				{
					sources_list.push_back(i);
				}
				N_source = mesh.vertices().size();
			}
			else
			{
				unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
				std::default_random_engine generator(seed);
				std::uniform_int_distribution<int> distribution(0, mesh.vertices().size() - 1);
				//std::srand((unsigned)time(NULL));
				unsigned s_i = 0;
				for (unsigned i = 0; i < N_source; i++)
				{
					//s_i = rand() % (mesh.vertices().size());
					s_i = distribution(generator);
					sources_list.push_back(s_i);
				}
			}
			fin.close();
			std::ofstream fout(sourceFile);
			for (unsigned i = 0; i < N_source; i++)
			{
				fout << sources_list[i] << std::endl;
			}
			fout.close();
		}
		else
		{
			unsigned s_id;
			while (fin >> s_id)
			{
				sources_list.push_back(s_id);
			}
			fin.close();
		}
		

	}
	else // single source
	{
		sources_list.clear();
		sources_list.push_back(source_id);
	}
	
	if (alg_id == 0) //-alg 0 -m -stxt (-o) (-t);   VTP
	{
		RunVTP();
	}
	else if (alg_id == 1) //-alg 1 -m -stxt -np -k (-o) (-t) (-judge);  PVTP
	{
		RunParallel_VTP();
	}
	else if (alg_id == 2) //-alg 2 -m -stxt -l (-o) (-t) (-c) (-r);  AVTP
	{
		RunAppr_VTP();
	}
	else if (alg_id == 3) //-alg 3 -m -stxt -l -np -k (-o) (-t) (-c) (-r). PAVTP
	{
		RunParallel_Appr_VTP();
	}
}


int main(int argc, char **argv)
{
	if (!InputParameter(argc, argv))
	{
		return -1;
	}
	
	std::vector<double> points;	
	std::vector<unsigned> faces;
	std::vector<int> realIndex;
	int originalVertNum = 0;
	
	std::cout << inFile << std::endl;
	
	// Load Mesh
	bool success = geodesic::read_mesh_from_file(inFile.c_str(), points, faces, realIndex, originalVertNum);
	if(!success)
	{
		std::cout << "something is wrong with the input file" << std::endl;
		return -2;
	}
	
	std::cout << "Load Mesh Success..." << std::endl;
	
	// Build Mesh
	mesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges

	if (mesh.IsDegenerate())
	{
		//system("pause");
		remove(inFile.c_str());
		return 0;
	}

	std::cout << "Build Mesh Success..." << std::endl;
	
	ComputeTau();

	RunAlgorithmByID(alg_id);

	
	//system("pause");
	return 0;
}