#include "itensor/all.h"
#include "cxxopts.hpp"
#include <boost/format.hpp>

/* #include "rang.hpp" */

/* #include <cmath> */
/* #include <vector> */
/* #include <sstream> */
/* #include <algorithm> */
/* #include <iomanip> */

using namespace itensor;
using namespace std;

/**************************************************************************//**
 *  Get num numbers between start and end (inclusive) that are equally spaced.
 * @param start The starting point of the range
 * @param end The end point of the range
 * @param num The number of points in the range
 * @return a vector of values with linear spacing
******************************************************************************/
template<typename T = double>
std::vector<T> linspace(T start, T end, int num) 
{
    T step = (end - start)/(num-1);
    T curValue = start;

    std::vector<T> vals;
    for (int n = 0; n < num; ++n) {
        vals.push_back(curValue);
        curValue += step;
    }
    return vals;
}

/**************************************************************************//**
 *  Get num numbers between 10^start and 10^end (inclusive) that appear 
 *  linearly spaced on a logarithmic scale.
 *
 * @param start The starting point of the range
 * @param end The end point of the range
 * @param num The number of points in the range
 * @return a vector of values with logarithmic spacing
******************************************************************************/
template<typename T = double>
std::vector<T> logspace(T start, T end, int num, T base = 10.0) 
{
    T step = (end - start)/(num-1);
    T curValue = start;

    std::vector<T> vals;
    for (int n = 0; n < num; ++n) {
        vals.push_back(pow(base, curValue));
        curValue += step;
    }
    return vals;
}

/**************************************************************************//**
 * MAIN PROGRAM
******************************************************************************/
int main(int argc, char* argv[]) {

    try {

        /* Store the command line options in a vector */
        vector<string> cmdl_args(argv,argv+argc); 

        cxxopts::Options options(argv[0], "Accessible Entanglement via DMRG");
        options
            .positional_help("[optional args]")
            .show_positional_help();

        options
            .allow_unrecognised_options()
            .add_options()
            ("pbc", "periodic boundary conditions")
            ("random", "use a random initial state")
            ("log_scale", "use a logarithmic spacing for the interaction strength")
            ("n,num_interaction", "the number of interaction points", cxxopts::value<int>(),"n")
            ("N, number_particles", "number of particles", cxxopts::value<int>(),"N")
            ("V, interaction_stength", "interaction strength", cxxopts::value<double>(),"V")
            ("i,Vi", "initial interaction strength", cxxopts::value<double>(),"Vi")
            ("f,Vf", "final interaction strength", cxxopts::value<double>(),"Vf")
            ("dV", "interaction strength step", cxxopts::value<double>(),"dV")
            ("help", "Print help");

        auto params = options.parse(argc, argv);

        if (params.count("help"))
        {
            cout << options.help({""}) << endl;
            exit(0);
        }
        if (!params.count("N"))
        {
            cout << options.help({""}) << endl;
            cout << "\nNeed to specify the number of particles N" << endl;
            exit(0);
        }
        if (!params.count("V")) {
            if (!(params.count("Vi") && params.count("Vf"))){
                cout << options.help({""}) << endl;
                cout << "\nNeed to specify the initial and final interaction strength Vi,Vf" << endl;
                exit(0);
            }
            /* else { */
            /*     cout << options.help({""}) << endl; */
            /*     cout << "\nNeed to specify at least one interaction strength V" << endl; */
            /*     exit(0); */
            /* } */
        }

        if ( (params.count("Vi") && params.count("Vf")) && 
             (params["Vi"].as<double>() > params["Vf"].as<double>()) ) {
            cout << options.help({""}) << endl;
            cout << "\n Vi must be smaller than Vf" << endl;
            exit(0);
        }

        if (!(params.count("dV") || params.count("num_interaction"))) {
            cout << options.help({""}) << endl;
            cout << "\nNeed to specify at least one of dV or num_interaction" << endl;
            exit(0);
        }
        if (params.count("dV") && params.count("num_interaction")) {
            cout << options.help({""}) << endl;
            cout << "\nCan't specify both dV and num_interaction" << endl;
            exit(0);
        }
        if (params.count("dV") && (params["dV"].as<double>() <= 0.0)) {
            cout << options.help({""}) << endl;
            cout << "\n dV must be positive" << endl;
            exit(0);
        }
        if (params.count("num_interaction") && (params["num_interaction"].as<int>() <= 0 )) {
            cout << options.help({""}) << endl;
            cout << "\n num_interaction must be a positive integer" << endl;
            exit(0);
        }

        double V,Vi,Vf;
        double dV,num_interaction;

        /* Collect and assign all options */
        if (params.count("V")){
            V = params["V"].as<double>();
            Vi = params["V"].as<double>();
            Vf = params["V"].as<double>();
        }
        else{
            V = params["Vi"].as<double>();
            Vi = params["Vi"].as<double>();
            Vf = params["Vf"].as<double>();
        }
        if (params.count("dV")) {
            dV = params["dV"].as<double>();
            num_interaction = floor((Vf-Vi)/dV) + 1;
        }
        else
            num_interaction = params["num_interaction"].as<int>();

        auto& N = params["N"].as<int>();
        auto L = 2*N;

        /* Now we create the list of interaction values */
        std::vector<double> interaction;
        if (params.count("log_scale")) {
            if (Vi > 0.0) 
                interaction = logspace(log10(Vi),log10(Vf),num_interaction);
            else {
                interaction = logspace(log10(abs(Vi)),log10(abs(Vf)),num_interaction);
                for(auto& cV : interaction)
                    cV *= -1;
            }
        
        }
        else 
            interaction = linspace(Vi,Vf,num_interaction);

        /* Do we include the Lth site?*/
        int last_site = L-1 + params.count("pbc");

        /* Create the output string */
        string out_file_name(str(boost::format("dmrg_acc_EE_%04d_%04d_%+11.5e_%+11.5e.dat") 
                % L % N % Vi % Vf));

        ofstream out_file;     // The output file
        out_file.open(out_file_name);

        /* Create the header of the file */
        out_file << "# ";
            for(string& carg: cmdl_args) 
                out_file << carg << " ";
        out_file << endl;

        /* string vN_label = str(boost::format("S₁(ℓ=%d)")%N); */
        string vN_label = str(boost::format("S1(l=%d)")%N);
        string vN_acc_label = str(boost::format("S1acc(l=%d)")%N);
        out_file << boost::format("#%17s%18s%18s%18s\n") % "V (t)" % "energy (t)" % vN_label % vN_acc_label;

        /* Create the sites basis of spinless fermions */
        int seed=7654;
        auto sites = Spinless(L,{"ConserveNf",true});

        //
        // Create the initial state
        // 
        auto psi = IQMPS(sites);
        char Occupation[2][5] = {"Emp", "Occ"}; 
        auto terms_element=InitState(sites);

        if(V >= 0.0)
        {
            for(int j = 0; j < 2; ++j)
            {
                int j_next=(j+1)%2;
                for(int i = 1; i <= L; ++i)
                {
                    if(i%2 == 0) terms_element.set(i,Occupation[j]);
                    else         terms_element.set(i,Occupation[j_next]);
                }
                if (j==0) psi= IQMPS(terms_element);
                else      psi=sum(psi,IQMPS(terms_element));   

            }
        }
        else if(V<0.0)
        {  
            for(int j = 1; j <= L; ++j)
            {
                int x=0;
                int y=1;
                int j_R=(j+L/2-2)%L+1;
                int j_L=j;
                if (j_R<j_L)
                {
                    j_L=j_R+1;
                    j_R=j-1;
                    x=1;
                    y=0;
                }
                for(int i = 1; i <= L; ++i)
                {
                    if((i>=j_L) && (i<=j_R)) terms_element.set(i,Occupation[y]);
                    else         terms_element.set(i,Occupation[x]);
                }
                if (j==1) psi= IQMPS(terms_element);
                else      psi=sum(psi,IQMPS(terms_element));   

            }
        }

        /* This code gives a seg fault! */
        /* if (params.count("random")) */
        /* { */
        /*     seedRNG(seed); */
        /*     psi = IQMPS(sites); */
        /*     psi.orthogonalize(); */
        /* } */

        for (const auto &cV : interaction) {

            /* Create the Hamiltonian */
            AutoMPO ampo(sites);
            for(int j = 1; j < L;++j)
            {
                ampo += -1.0,"Cdag",j,"C",j+1;
                ampo += -1.0,"Cdag",j+1,"C",j;
                ampo += cV,"N",j+1,"N",j;
            }

            if (params.count("pbc")) {
                auto factor = 1.0;
                if (L/2 % 2 == 0)
                    factor *= -1;
                ampo += -1.0*factor,"Cdag",1,"C",L;
                ampo += -1.0*factor,"Cdag",L,"C",1;
                ampo += cV,"N",1,"N",L;
            }
            auto H = IQMPO(ampo);

            auto sweeps = Sweeps(5);
            sweeps.maxm() = 10,20,100,100,200;
            sweeps.cutoff() = 1E-10;
            sweeps.niter() = 2;
            sweeps.noise() = 1E-7,1E-8,0.0;

            /* Perform the DMRG calculation */
            normalize(psi);
            auto energy = dmrg(psi,H,sweeps,"Quiet");

            /* Compute the vN entanglement entropy */
            auto i = L/2;
            psi.position(i);
            auto wf = psi.A(i)*psi.A(i+1);
            auto U = psi.A(i);
            IQTensor S,VV;
            auto spectrum = svd(wf,U,S,VV);

            Real SvN = 0.;
            for(auto p : spectrum.eigs())
            {
                if(p > 1E-12) SvN += -p*log(p);
            }

            auto firstofS=S.index(1);
            int Nn=firstofS.nblock();
            vector<Real> Pn(Nn);
            vector<Real> Pna(Nn);
            int counter=0;

            for (int i=0; i<Nn; i++)
            {
                for (int j=1; j<= firstofS.index(i+1).m(); j++)
                {
                    counter+=1;
                    Pn[i]+= pow(S.real(counter,counter),2.0);
                    Pna[i]+= pow(S.real(counter,counter),4.0);
                }
            }
            Real PnSum=0;
            Real PnaSum=0;
            for (int i=0; i<Nn; i++)
            {
                PnSum+=Pn[i];
                PnaSum +=Pna[i];
            }
            for (int i=0; i<Nn; i++)
            {
                Pn[i]/= PnSum;
                Pna[i]/= PnaSum;
            }
            Real PnvN = 0.;

            for(auto pn : Pn)
            {
                if(pn > 1E-12) PnvN += -pn*log(pn);
            }
            Real PnavN = 0.;

            for(auto pna : Pna)
            {
                if(pna > 1E-12) PnavN += -pna*log(pna);
            }

            /* printfln("%18.8e%18.8e%18.8e\n",energy,SvN,(SvN-PnvN)); */
            out_file << boost::format("%18.8e%18.8e%18.8e%18.8e\n") % cV % energy % SvN % (SvN-PnvN);
            out_file.flush();

        }

        /* Close the output file */
        if (out_file.is_open())
            out_file.close();

    } 
    catch (const cxxopts::OptionException& e)
    {
        cout << "error parsing options: " << e.what() << endl;
        exit(1);
    }

    return 0;
}


