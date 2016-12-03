// File:    mmsp2xyz.cpp
// Purpose: reads MMSP grid containing sparse floats, tracks specified grain
// Output:  CSV file specifying XYZ+phase
// Depends: MMSP, zlib

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <zlib.h>
#include <sstream>
#include <cmath>
#include <map>
#include "MMSP.hpp"
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <stdlib.h>

using namespace MMSP;

int main(int argc, char* argv[]) {
	if ( argc != 2 && argc != 3 && argc != 6) {
		std::cout << "Usage: " << argv[0] << " xxx.dat [idToPixelNumberMap->txt] [x_low_bound] [x_high_bound [dim_of_grad]]\n";
		return ( 1 );
	}
	
	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[1] << ".\n\n";
		exit(-1);
	}

	// read data type
	std::string type;
	getline(input, type, '\n');

	int gradDir = 0;
	if (argc == 6) {
	  gradDir = atoi(argv[5]);
	}
	std::cout << "argc:" << argc << "  gradDir: " << gradDir << std::endl;

  bool unsigned_long_type = (type.find("unsigned long") != std::string::npos);
  if (unsigned_long_type == false) {
    std::cerr << "Only output with unsigned long type is supported, your data type is " << type.substr(type.find(":") + 1) << ".\n\n";
    exit(-1);
  }    


	// read grid dimension
  int dim;
  input >> dim;
  #ifdef DEBUG
  std::cout<<"Grid is "<<dim<<"-dimensional."<<std::endl;
  #endif

  // read number of fields
  int fields;
  input >> fields;
  #ifdef DEBUG
  std::cout<<"Grid has "<<fields<<" fields."<<std::endl;
  #endif

  // read grid sizes
  int g0[3] = {0, 0, 0};
  int g1[3] = {0, 0, 0};
  for (int i = 0; i < dim; i++)
	  input >> g0[i] >> g1[i];
  #ifdef DEBUG
  std::cout<<"Grid edge is " <<g1[0] - g0[0]<<std::endl;
  #endif

  // read cell spacing
  float dx[3] = {1.0, 1.0, 1.0};
  for (int i = 0; i < dim; i++)
	  input >> dx[i];
  #ifdef DEBUG
  std::cout << "Grid spacing is " << dx[0] << std::endl;
  #endif

  // ignore trailing endlines
  input.ignore(10, '\n');

  // determine byte order
  std::string byte_order;
  if (0x01 & static_cast<int>(1)) {
    byte_order = "LittleEndian";
  }
  else {
    byte_order = "BigEndian";
  }
  #ifdef DEBUG
  std::cout << "Endian is "<< byte_order <<std::endl;
  #endif

  // read number of blocks
  int blocks;
  input.read(reinterpret_cast<char*>(&blocks), sizeof(blocks));
  std::unordered_set<unsigned long> *boundaryGrainId = new std::unordered_set<unsigned long>;
  std::unordered_map<unsigned long, unsigned long> *idToPixelNumberMap = new std::unordered_map<unsigned long, unsigned long>; // only consider interior grains 

  #ifdef DEBUG
  std::cout << "container generated " <<std::endl;
  #endif

  unsigned long length = input.tellg();

  // we scan through the grid twice. 
  // first scan to find boundary grain id
  for (int i = 0; i < blocks; i++) {
	  // read block limits
	  int lmin[3] = {0, 0, 0};
	  int lmax[3] = {0, 0, 0};
    int gmin[3] = {0, 0, 0};
    int gmax[3] = {0, 0, 0};

    for (int j = 0; j < dim; j++) {
      gmin[j] = g0[j];
      gmax[j] = g1[j] - 1;
    }

    if (argc == 6) {
      gmin[gradDir] = atoi(argv[3]);
      gmax[gradDir] = atoi(argv[4]);
      #ifdef DEBUG
        std::cout << "dim: " << gradDir << "   min: " << gmin[gradDir] << " max: " << gmax[gradDir] << std::endl;
      #endif
    }


	  for (int j = 0; j < dim; j++) {
		  input.read(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
		  input.read(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
	  }
	  //#ifdef DEBUG
	  //std::cout<<"  Block edge is "<<lmax[0] - lmin[0]<<std::endl;
	  //#endif
	  int blo[dim];
    int bhi[dim];
    // read boundary conditions
    for (int j = 0; j < dim; j++) {
      input.read(reinterpret_cast<char*>(&blo[j]), sizeof(blo[j]));
      input.read(reinterpret_cast<char*>(&bhi[j]), sizeof(bhi[j]));
    }

	  // read grid data
	  unsigned long size, rawSize;
	  input.read(reinterpret_cast<char*>(&rawSize), sizeof(rawSize)); // read raw size
	  input.read(reinterpret_cast<char*>(&size), sizeof(size)); // read compressed size
	  char* compressed_buffer = new char[size];
	  input.read(compressed_buffer, size);
	  //#ifdef DEBUG
	  //std::cout<<"  Read "<<size<<" B, compressed data."<<std::endl;
	  //#endif
	  char* buffer;
		if (size!=rawSize) {
		  // Decompress data
		  buffer = new char[rawSize];
		  int status;
		  status = uncompress(reinterpret_cast<unsigned char*>(buffer), &rawSize, reinterpret_cast<unsigned char*>(compressed_buffer), size);
		  switch( status ) {
		  case Z_OK:
			  break;
		  case Z_MEM_ERROR:
			  std::cerr << "Uncompress: out of memory." << std::endl;
			  exit(1);
			  break;
		  case Z_BUF_ERROR:
			  std::cerr << "Uncompress: output buffer wasn't large enough." << std::endl;
			  exit(1);
			  break;
		  }
		  delete [] compressed_buffer;
	  } else {
		  buffer=compressed_buffer;
		  compressed_buffer=NULL;
	  }


    if (dim == 2) {
			MMSP::grid<2,MMSP::scalar<unsigned long> > GRID(fields, lmin, lmax);
		  GRID.from_buffer(buffer);
		  for (int k = 0; k < (lmax[1]-lmin[1]); k++){
        for (int l = 0; l < (lmax[0]-lmin[0]); l++){
          vector<int> x (2,0);
          x[0] = lmin[0] + l;
          x[1] = lmin[1] + k;
          if (argc == 6 && (x[gradDir] < gmin[gradDir] || x[gradDir] > gmax[gradDir])) {
            continue;
          }
          if (x[0] == gmin[0] || x[0] == gmax[0] || x[1] == gmin[1] || x[1] == gmax[1]) {
            if (boundaryGrainId->find(GRID(x)) == boundaryGrainId->end()) {
              boundaryGrainId->insert(GRID(x));
            }
          }
        } // for l
		  } // for k 	
		} else if (dim == 3) {
  	  MMSP::grid<3,MMSP::scalar<unsigned long> > GRID(fields, lmin, lmax);
		  GRID.from_buffer(buffer);
      for (int l = 0; l < (lmax[0]-lmin[0]); l++){
		    for (int k = 0; k < (lmax[1]-lmin[1]); k++){
          for (int m = 0; m < (lmax[2]-lmin[2]); m++){
            vector<int> x (3,0);
            x[0] = lmin[0] + l;
            x[1] = lmin[1] + k;
            x[2] = lmin[2] + m;
            if (argc == 6 && (x[gradDir] < gmin[gradDir] || x[gradDir] > gmax[gradDir])) {
              continue;
            }
            if (x[0] == gmin[0] || x[0] == gmax[0] || x[1] == gmin[1] || x[1] == gmax[1] || x[2] == gmin[2] || x[2] == gmax[2]) {
              if (boundaryGrainId->find(GRID(x)) == boundaryGrainId->end()) {
                boundaryGrainId->insert(GRID(x));
              }
            }
          } // for m
        } // for k
      }// for l 
		}// if dim == 3
	  // clean up
	  delete [] buffer; 
    buffer = NULL;
	}// for int i //loop over blocks

  #ifdef DEBUG
  std::cout << "first scan finished" <<std::endl;
  #endif

  input.seekg (length, input.beg);

  // scan the second time for id to #pixel map
  for (int i = 0; i < blocks; i++) {
    // read block limits
    int lmin[3] = {0, 0, 0};
    int lmax[3] = {0, 0, 0};

    int gmin[3] = {0, 0, 0};
    int gmax[3] = {0, 0, 0};

    for (int j = 0; j < dim; j++) {
      gmin[j] = g0[j];
      gmax[j] = g1[j] - 1;
    }

    if (argc == 6) {
      gmin[gradDir] = atoi(argv[3]);
      gmax[gradDir] = atoi(argv[4]);
      #ifdef DEBUG
        std::cout << "dim: " << gradDir << "   min: " << gmin[gradDir] << " max: " << gmax[gradDir] << std::endl;
      #endif
    }

    for (int j = 0; j < dim; j++) {
      input.read(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
      input.read(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
    }
    //#ifdef DEBUG
    //std::cout<<"  Block edge is "<<lmax[0] - lmin[0]<<std::endl;
    //#endif
    int blo[dim];
    int bhi[dim];
    // read boundary conditions
    for (int j = 0; j < dim; j++) {
      input.read(reinterpret_cast<char*>(&blo[j]), sizeof(blo[j]));
      input.read(reinterpret_cast<char*>(&bhi[j]), sizeof(bhi[j]));
    }

    // read grid data
    unsigned long size, rawSize;
    input.read(reinterpret_cast<char*>(&rawSize), sizeof(rawSize)); // read raw size
    input.read(reinterpret_cast<char*>(&size), sizeof(size)); // read compressed size
    char* compressed_buffer = new char[size];
    input.read(compressed_buffer, size);
    //#ifdef DEBUG
    //std::cout<<"  Read "<<size<<" B, compressed data."<<std::endl;
    //#endif
    char* buffer;
    if (size!=rawSize) {
      // Decompress data
      buffer = new char[rawSize];
      int status;
      status = uncompress(reinterpret_cast<unsigned char*>(buffer), &rawSize, reinterpret_cast<unsigned char*>(compressed_buffer), size);
      switch( status ) {
      case Z_OK:
        break;
      case Z_MEM_ERROR:
        std::cerr << "Uncompress: out of memory." << std::endl;
        exit(1);
        break;
      case Z_BUF_ERROR:
        std::cerr << "Uncompress: output buffer wasn't large enough." << std::endl;
        exit(1);
        break;
      }
      delete [] compressed_buffer;
    } else {
      buffer=compressed_buffer;
      compressed_buffer=NULL;
    }


    if (dim == 2) {
      MMSP::grid<2,MMSP::scalar<unsigned long> > GRID(fields, lmin, lmax);
      GRID.from_buffer(buffer);
      for (int k = 0; k < (lmax[1]-lmin[1]); k++){
        for (int l = 0; l < (lmax[0]-lmin[0]); l++){
          vector<int> x (2,0);
          x[0] = lmin[0] + l;
          x[1] = lmin[1] + k;
          if (argc == 6 && (x[gradDir] < gmin[gradDir] || x[gradDir] > gmax[gradDir])) {
            continue;
          }
          unsigned long id = GRID(x);
          if (boundaryGrainId->find(id) == boundaryGrainId->end()) {
            if (idToPixelNumberMap->find(id) == idToPixelNumberMap->end()) {
              idToPixelNumberMap->insert({id, (unsigned long)1});
            } else {
              // because [] is overloaded in scope of MMSP, so we use find() here;
              (idToPixelNumberMap->find(id)->second)++;
            }  
          }
        } // for l
      } // for k  
    } else if (dim == 3) {
      MMSP::grid<3,MMSP::scalar<unsigned long> > GRID(fields, lmin, lmax);
      GRID.from_buffer(buffer);
      for (int l = 0; l < (lmax[0]-lmin[0]); l++){
        for (int k = 0; k < (lmax[1]-lmin[1]); k++){
          for (int m = 0; m < (lmax[2]-lmin[2]); m++){
            vector<int> x (3,0);
            x[0] = lmin[0] + l;
            x[1] = lmin[1] + k;
            x[2] = lmin[2] + m;
            if (argc == 6 && (x[gradDir] < gmin[gradDir] || x[gradDir] > gmax[gradDir])) {
              continue;
            }
            unsigned long id = GRID(x);
            if (boundaryGrainId->find(id) == boundaryGrainId->end()) {
              if (idToPixelNumberMap->find(id) == idToPixelNumberMap->end()) {
                idToPixelNumberMap->insert({id, (unsigned long)1}); // map insert must be in pair format {}
              } else {
                // because [] is overloaded in scope of MMSP, so we use find() here;
                (idToPixelNumberMap->find(id)->second)++;
              }  
            }
          } // for m
        } // for k
      }// for l 
    }// if dim == 3
    // clean up
    delete [] buffer; 
    buffer = NULL;
  }// for int i //loop over blocks

  #ifdef DEBUG
  std::cout << "Second scan finished " <<std::endl;
  #endif

  //std::cout << "interior grain number: " << idToPixelNumberMap->size() << std::endl;
  //std::cout << "boundary grain number: " << boundaryGrainId->size() << std::endl;
  std::string str(argv[1]);
  std::string str1 = str.substr(str.find_first_of('.') + 1);
  std::cout << str1.substr(0, str1.find_first_of('.')) << "  " <<  idToPixelNumberMap->size() << "  " << boundaryGrainId->size() << std::endl; 

  if (argc == 3 || argc == 6) {
    std::ofstream outFile;
    std::string outputfile(argv[2]);
    outFile.open(outputfile);

    for (std::unordered_map<unsigned long, unsigned long>::iterator i = idToPixelNumberMap->begin(); i != idToPixelNumberMap->end(); i++) {
      outFile << i->first << " " << i->second << "\n";
    }

    outFile.close();
  }
  

  #ifdef DEBUG
  std::cout << "Output finished " <<std::endl;
  #endif

  delete boundaryGrainId;
  delete idToPixelNumberMap;
	return 0;
}
