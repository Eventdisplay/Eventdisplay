/*! \file mergeVBF.cpp
    \brief merge several vbf files into one

*/


#include "VBankFileReader.h"
#include "VBankFileWriter.h"
#include "VPacket.h"
#include "VArrayEvent.h"
#include "VDatum.h"

// include the simulation data structure
#include "VSimulationData.h"

// include the configuration mask utilities, which give us parseConfigMask()
#include "VConfigMaskUtil.h"


#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

using namespace VConfigMaskUtil;


void usage(char *prog){
  cout<<"Usage: "<<prog<<" [listOfFiles] [output.vbf] [newRunNumber]"<<endl;
  exit(-1);
}

int main(int argc, char **argv){
  
  if(argc!=4)
    usage(argv[0]);

  vector<string> fileNames;
  
  ifstream infile;
  infile.open(argv[1]);
  if(!infile.is_open()){
    cerr<<"Couldn't open "<<argv[1]<<endl;
    exit(-1);
  }

  int newRunNumber=0;
  sscanf(argv[3], "%d", &newRunNumber);

  while(true){
    char buffer[1000];
    infile.getline(buffer, 1000);
    if(infile.eof())
      break;
    fileNames.push_back(buffer);
  }

  cout << "Collected " << fileNames.size() << " files" << endl;

  try{  

    VBankFileWriter writer(argv[2], newRunNumber, parseConfigMask("0,1,2,3"));

    VPacket *packet=NULL;
    VArrayTrigger *trigger=NULL;
    VArrayEvent *arrayEvent=NULL;
    VSimulationData *sim=NULL;
    VSimulationHeader *header=NULL;
    VEvent *telEvent = NULL;
  
    bool wroteSimHeader=false;

    int globalEventCount=0;

    bool writePacket=true;

    for (unsigned fileIndex=0; fileIndex<fileNames.size(); fileIndex++)
      {
	cout << "reading file " << fileIndex << ": " << fileNames[fileIndex] << endl;

	try{

	  //VBankFileReader reader(fileNames.at(fileIndex).c_str(), false, false );
	  VBankFileReader reader(fileNames.at(fileIndex).c_str() );
	  //reader.generateIndexAndChecksum();
	  int numPackets = reader.numPackets();
	  cout<<"\t Packets: "<<numPackets<<endl;
    
      
	  for (int pack=0; pack<numPackets; pack++){
	    packet = reader.readPacket(pack);
	    if(packet){
	        
	      writePacket=true;
	        
	      if(packet->hasArrayEvent()){
		arrayEvent = packet->getArrayEvent();
		if(arrayEvent){
		  trigger = arrayEvent->getTrigger();
		  if(trigger){
		    trigger->setRunNumber(newRunNumber);
		    trigger->setEventNumber(globalEventCount);
		  }
		  for (unsigned i=0; i<arrayEvent->getNumEvents(); i++){
		    telEvent = arrayEvent->getEvent(i);
		    if(telEvent){
		      telEvent->setEventNumber(globalEventCount);
		    }
		  }
		}
	      }
	        
	      if(packet->hasSimulationData()){
		sim = packet->getSimulationData();
		if(sim){
		  sim->fRunNumber = newRunNumber;
		  sim->fEventNumber=  globalEventCount;
		}
	      }
	        
	      if(packet->hasSimulationHeader()){
		header = packet->getSimulationHeader();
		if(header){
		  header->fRunNumber = newRunNumber;
		}
		    
		if(wroteSimHeader)
		  writePacket=false;
		    
		wroteSimHeader=true;
		    
	      }
	        
	      if(writePacket){
		writer.writePacket(packet);
		globalEventCount++;
	      }
	        
	      delete packet;
	    }
    
	  }
	}
	catch (const std::exception&ex)
	  {
	    cerr <<"For file "<<fileNames.at(fileIndex)<<endl;
	    cout << ex.what() << endl;
	    //  cerr <<ex<<endl;
	  }
      } 
    writer.finish();
  }
  catch (const std::exception &ex)
    {
      //      cerr<<ex<<endl;
      exit(EXIT_FAILURE);
    }
  /*  catch (const std::exception &ex)
    {
    //      cerr<<ex<<endl;
      exit(EXIT_FAILURE);
    }
  catch (const std::exception& e)
    {
    //      cerr<<e.what()<<endl;
      exit(EXIT_FAILURE);
      } */
  catch(...)
    {
      //      cerr<<"Unknown exception caught! Program crashing and burning"<<endl;
      exit(EXIT_FAILURE);
    } 
  
}
