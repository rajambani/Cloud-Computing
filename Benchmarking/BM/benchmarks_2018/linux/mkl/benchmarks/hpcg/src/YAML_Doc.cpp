/*******************************************************************************
* Copyright 2014-2017 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/


//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include "YAML_Doc.hpp"
using namespace std;

/*!
  Sets the application name and version which will become part of the YAML doc.

  @param[in] miniApp_Name application name
  @param[in] miniApp_Version application name
  @param[in] destination_Directory destination directory for the YAML document
  @param[in] destination_FileName file name for the YAML document
*/
YAML_Doc::YAML_Doc(const std::string & miniApp_Name, const std::string & miniApp_Version, const std::string & destination_Directory, const std::string & destination_FileName) {
  miniAppName = miniApp_Name;
  miniAppVersion = miniApp_Version;
  destinationDirectory = destination_Directory;
  destinationFileName = destination_FileName;
}

//inherits the destructor from YAML_Element
YAML_Doc::~YAML_Doc(void) {
}

/*!
  Generates YAML from the elements of the document and saves it to a file.

  @return returns the complete YAML document as a string
*/
string YAML_Doc::generateYAML() {
  string yaml;

  yaml =  yaml + miniAppName + " version: " + miniAppVersion + "\n";

  for (size_t i=0; i<children.size(); i++) {
    yaml = yaml + children[i]->printYAML("");
  }

  time_t rawtime;
  tm * ptm;
  time ( &rawtime );
  ptm = localtime(&rawtime);
  char sdate[25];
  //use tm_mon+1 because tm_mon is 0 .. 11 instead of 1 .. 12
  sprintf (sdate,"%04d.%02d.%02d.%02d.%02d.%02d",ptm->tm_year + 1900, ptm->tm_mon+1,
      ptm->tm_mday, ptm->tm_hour, ptm->tm_min,ptm->tm_sec);

  string filename;
  if (destinationFileName=="") {
    filename = miniAppName + "-" + miniAppVersion + "_";
    filename = filename + string(sdate) + ".yaml";
  }
  else {
    filename = destinationFileName;
    filename = filename + ".yaml";
  }
  

  if (destinationDirectory!="" && destinationDirectory!=".") {
    string mkdir_cmd = "mkdir " + destinationDirectory;
    system(mkdir_cmd.c_str());
    filename = destinationDirectory + "/" + destinationFileName;
  } else
    filename = "./" + filename;

  ofstream myfile;
  myfile.open(filename.c_str());
  myfile << yaml;
  myfile.close();
  return yaml;
}
