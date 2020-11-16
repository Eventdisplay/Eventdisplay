/*! \file splitVBF.cpp
    \brief split one vbf file into several ones

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
  cout << "Usage: " << prog << " [input.vbf] [numberOfFiles] [newRunNumber]" << endl;
  exit(-1);
}

int main(int argc, char **argv){
  
  if(argc != 4)
    usage(argv[0]);

  int numberOfFiles = 0;
  sscanf(argv[2], "%d", &numberOfFiles);

  cout << "NumberOfFiles: " << numberOfFiles << endl;

  int newRunNumber = 1;
  sscanf(argv[3], "%d", &newRunNumber);

  cout << "newRunNumber: " << newRunNumber << endl;

  cout << "Reading input file " << argv[1] << endl;
  
  VBankFileReader reader(argv[1]);
  const int numPackets = reader.numPackets();
  cout << "Packets: " << numPackets << endl;

  const int numPacketsPerFile = int(numPackets / numberOfFiles);
  cout << "\nPackets/file: " << numPacketsPerFile << endl;

  // Extraction the .vbf from the input filename
  string inputFileName(argv[1]);
  const size_t dot_index = inputFileName.find_last_of("."); 
  const string rawInputName = inputFileName.substr(0, dot_index);

  cout << "Raw filenames: " << rawInputName << endl;

  VPacket *packet = NULL;
  VPacket *packetHeader = NULL;
  VArrayTrigger *trigger = NULL;
  VArrayEvent *arrayEvent = NULL;
  VSimulationData *sim = NULL;
  VSimulationHeader *header = NULL;
  VEvent *telEvent = NULL;

  try{
      // Searching for sim header and storing it
      int ipacketHeader = -1;
      for (int ipacket = 0; ipacket < numPackets; ipacket++) {
        packetHeader = reader.readPacket(ipacket);
        if(packetHeader && packetHeader->hasSimulationHeader()){
          cout << "Simulation header found at packet " << ipacket << endl;
          header = packetHeader->getSimulationHeader();          
          header->fRunNumber = newRunNumber;
          ipacketHeader = ipacket;
          break;        
        } // if packet
      } // ipacket

      for (int ifile = 0; ifile < numberOfFiles; ifile++) {
        
        stringstream fileNameTmp;
        fileNameTmp << rawInputName << "_" << ifile + 1 << ".vbf";
        string fileName = fileNameTmp.str();

        cout << "ifile = " << ifile << endl;
        cout << "FileName: " << fileName << endl;

        VBankFileWriter writer(fileName, newRunNumber, parseConfigMask("0,1,2,3"));
      
        int globalEventCount = 1;
        bool writePacket = true;

        const int startPacket = ifile * numPacketsPerFile;
        const int endPacket = (ifile + 1) * numPacketsPerFile;

        if (header) {
          cout << "Writing sim header" << endl;
          writer.writePacket(packetHeader);
        }

        for (int ipacket = startPacket; ipacket < endPacket; ipacket++) {

          if (ipacket == ipacketHeader)
            continue;

          writePacket = false;

          packet = reader.readPacket(ipacket);

          // cout << "Starting ipacket: " << ipacket << endl;

          if(packet){
            writePacket = true;
                
            if(packet->hasArrayEvent()){
              arrayEvent = packet->getArrayEvent();
            
              // cout << "Packet " << ipacket << " has arrayEvent" << endl;

              if(arrayEvent){
                trigger = arrayEvent->getTrigger();
              
                if(trigger){
                  trigger->setRunNumber(newRunNumber);
                  trigger->setEventNumber(globalEventCount);
                }
              
                for (unsigned ievent = 0; ievent < arrayEvent->getNumEvents(); ievent++){
                  telEvent = arrayEvent->getEvent(ievent);
                      
                  if(telEvent){
                    telEvent->setEventNumber(globalEventCount);
                  }
                } // ievent
              } // arrayEvent
            } // packet
                
            if(packet->hasSimulationData()){
              sim = packet->getSimulationData();
              
              if(sim){
                sim->fRunNumber = newRunNumber;
                sim->fEventNumber=  globalEventCount;
              }
            } //hasSimulationData
                            
            if(writePacket){
              writer.writePacket(packet);
              globalEventCount++;
            }
        
          } // if packet
            
          delete packet;
        } // ipacket

        writer.finish();
      } // ifile 

      delete packetHeader;
  } // try
  catch(...)
      {
        exit(EXIT_FAILURE);
      } 
    
}
