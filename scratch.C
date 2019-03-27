
    else if (settings1->NOISE == 2) {    // NOISE == 2 : use Rayleigh dist. fit from A2/A3 data

        Vfft_noise_after.clear();  // remove previous Vfft_noise values
        Vfft_noise_before.clear();  // remove previous Vfft_noise values
        //V_noise_timedomain.clear(); // remove previous V_noise_timedomain values


        double V_tmp; // copy original flat H_n [V] value
        double current_amplitude, current_phase;

        GetNoisePhase(settings1); // get random phase for noise
	std::ifstream infile("sigmavsfreq_ch0.txt");
	double  a, b, c, d;
	vector <double> freqs;
	vector <double> V_rayleigh;
	freqs.clear();
	while (infile >> a >> b >> c >> d)
	  {
	    TF1 *f1 = new TF1("f1","([0]*x)/([1]*[1])*exp(-[0]*[0]*x*x/(2.*[1]*[1]))");  
	    f1->SetParameters(b,c);
	    double r = f1->GetRandom();
	    freqs.push_back(a);
	    V_rayleigh.push_back(r);

	    for (int k=0; k<settings1->DATA_BIN_SIZE/2; k++) {


	      current_phase = noise_phase[k];
	      
	      
	      //   cout<<"Random number from Rayleigh dist "<<r<<endl;
	    }
	    //  TCanvas *c1 = new TCanvas("","",850,850);
	    
	    vnoise[2 * k] = V_rayleigh[k] * cos(noise_phase[k]);
	    vnoise[2 * k + 1] = V_rayleigh[k] * sin(noise_phase[k]);


	    Vfft_noise_after.push_back( vnoise[2*k] );
	    Vfft_noise_after.push_back( vnoise[2*k+1] );

	    // inverse FFT normalization factor!
	    vnoise[2 * k] *= 2./((double)settings1->DATA_BIN_SIZE);
	    vnoise[2 * k + 1] *= 2./((double)settings1->DATA_BIN_SIZE);
	    

        }



        // now vnoise is time domain waveform
        Tools::realft( vnoise, -1, settings1->DATA_BIN_SIZE);

     
    }
