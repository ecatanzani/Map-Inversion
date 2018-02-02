
#include "MyHead.h"

string output_path_creator(const Int_t out_choose) {
  // out_choose == 0 means we are creating the path for a log file
  // out_choose == 0 means we are creating the path for a ROOT output file

  string output;
  
  switch(out_choose) {

  case 0:
    output = output_log;
    output += string_sbi_subsample;
    output += "_";
    output+=to_string((long long)time_stamp);
    output+=".txt";
    cout<<"\nWritten log file: -> \t "<<output<<endl;
    break;

  case 1:
    output = output_root;
    output += string_sbi_subsample;
    output += "_";
    output+=to_string((long long)time_stamp);
    output+="_maps_result.root";
    cout<<"\nWritten ROOT file: -> \t "<<output<<endl;
    break;

  }
  
  return output;
}

void log_file_init(ofstream &out_file) {
  out_file << "********************* Automatic Log File Generator *******************"<<endl<<endl;
  
  out_file << "////////////////////////// Simulation Parameters //////////////////////////"<<endl<<endl;
  out_file << "Simulation timestamp: "<<time_stamp<<endl;
  out_file << "Simulation TRandom3 seed: "<<random_seed<<endl;
  out_file << "Number of readed SBI files: "<<number_SBI_files<<endl;
  out_file << "Number of generated events per sky portion: "<<sky_events<<endl<<endl;
 
  out_file << "*\n*\n*\n"<<endl;
}
