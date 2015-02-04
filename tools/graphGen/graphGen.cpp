/**
 * \file graphGen.cpp
 * \brief Generator for sample graph data used in system tests. See also README.
 * \author Gregor Matura
 * \version 0.1.0
 * \date 2013-05-14
 */
/* \copyright Copyright 2013 German Aerospace Center (http://www.DLR.de)\n\n
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iterator>

#define VEC_LEN_START "["
#define VEC_LEN_END "]"
#define VEC_START "{"
#define VEC_END "}"
#define SEPARATOR ", "

#define STRUCT_NAME "exGraphData"
#define ARRAY_TYPE_DEF "int"
#define ARRAY_TYPE_DECLARE "static int"

#define WGT_BASE 1000
#define MOD_WGTCMP 100
#define MOD_WGTVIS 200
#define FAC_PROC 10
#define FAC_PART 1
#define FAC_SYM FAC_PART

#define DEF_WIDTH 4

using std::cout;
using std::endl;
using std::stringstream;
using std::string;
using std::vector;

class rect {
  public:
    int w;
    int h;
  
    rect(int w, int h) { this->w = w; this->h = h; };
    
    int getLinIdxAt(int i, int j) { return j*w+i; };
};

template <typename T>
string vecToStr(vector<T> vec) {
  stringstream sstr("");
  sstr << VEC_START;
  for (typename vector<T>::iterator it=vec.begin(); it!=vec.end(); ++it) {
    sstr << (it!=vec.begin() ? SEPARATOR : "") << std::setw(DEF_WIDTH) << *it;
  }
  sstr << VEC_END;
  return sstr.str();
}

template <typename T>
string vecInfoSizeOnly(vector<T> vec, string description1="", string description2="") {
  stringstream sstr;
  sstr << description1 << description2;
  sstr << VEC_LEN_START << vec.size() << VEC_LEN_END;
  return sstr.str();
}

template <typename T>
string vecInfo(vector<T> vec, string description1="", string description2="") {
  stringstream sstr;
  sstr << vecInfoSizeOnly<T>(vec, description1, description2);
  sstr << " = ";
  sstr << vecToStr<T>(vec);
  return sstr.str();
}

template<typename T>
T strToNbr(string str) {
  stringstream sstr(str);
  T ret;
  sstr >> ret;
  return ret;
}

void genGlobal(int w, int h, vector<int>& vertglb, vector<int>& edgeglb) {
  rect mRect(w,h);

  vertglb.push_back(0);

  for (int j=0; j<mRect.h; ++j) {
    for (int i=0; i<mRect.w; ++i) {
      if (j>0) { edgeglb.push_back(mRect.getLinIdxAt(i, j-1)); }
      if (i>0) { edgeglb.push_back(mRect.getLinIdxAt(i-1, j)); }
      if (i<mRect.w-1) { edgeglb.push_back(mRect.getLinIdxAt(i+1, j)); }
      if (j<mRect.h-1) { edgeglb.push_back(mRect.getLinIdxAt(i, j+1)); }
      vertglb.push_back(edgeglb.size());
    }
  }
  
}

void getPart(int numParts, int curPart, const vector<int>& vertglb, const vector<int>& edgeglb, vector<int>& vertloc, vector<int>& edgeloc) {
    int chStart = curPart * vertglb.size() / numParts;
    int chEnd = std::min(vertglb.size(), (curPart+1)*vertglb.size()/numParts+1);

    vertloc.insert(vertloc.begin(), vertglb.begin()+chStart, vertglb.begin()+chEnd);

    edgeloc.insert(edgeloc.begin(), edgeglb.begin()+vertloc[0], edgeglb.begin()+vertloc[vertloc.size()-1]);

    for (vector<int>::reverse_iterator rit=vertloc.rbegin(); rit!=vertloc.rend(); ++rit) {
      *rit -= vertloc[0];
    }
}

template <typename T>
vector<T> concatVecVec(vector<vector<T> > in) {
  vector<T> ret;
  for (typename vector<vector<T> >::iterator it=in.begin(); it!=in.end(); ++it) {
    ret.insert(ret.end(), it->begin(), it->end());
  }
  return ret;
}


vector<vector<int> > genWgtsForParts(vector<vector<int> > in, int base = WGT_BASE) {
  vector<vector<int> > ret;
  for (vector<vector<int> >::iterator outerIt=in.begin(); outerIt!=in.end(); ++outerIt) {
    vector<int> curPart;
    for (vector<int>::iterator innerIt=outerIt->begin(); innerIt!=outerIt->end(); ++innerIt) {
      curPart.push_back(base + FAC_PROC*(outerIt - in.begin()) + FAC_PART*(innerIt - outerIt->begin()));
    }
    ret.push_back(curPart);
  }
  return ret;
}

vector<vector<int> > genConWgtsForParts(vector<vector<int> > in, int base = WGT_BASE) {
  vector<vector<int> > ret;
  for (vector<vector<int> >::iterator outerIt=in.begin(); outerIt!=in.end(); ++outerIt) {
    vector<int> curPart;
    for (vector<int>::iterator innerIt=outerIt->begin(); innerIt!=outerIt->end(); ++innerIt) {
      curPart.push_back(base);
    }
    ret.push_back(curPart);
  }
  return ret;
}

int findVertIndex(vector<int> const& in, int const needle) {
  std::vector<int>::const_iterator it = in.begin() + 1;
  while (it!=in.end()) {
    if (needle < *it) return it - in.begin() - 1;
    ++it;
  }
  return in.end() - in.begin() - 1;
}

/**
 * \brief Get vertex index of edge tail.
 *
 * Find edge index in appropriate (=eProc) part of partsVert.
 */
int getTailVertexIndex(vector<int> const& vtx, vector<vector<int> > const& partsVert, vector<vector<int> > const& partsEdge, int const eProc, int const eIdx) {
  int vLocIdx = findVertIndex(partsVert[eProc], eIdx);
  int vGlbOffset = vtx[eProc+1-1];
  return vGlbOffset + vLocIdx;
}

void findPairingEdgeIndices(vector<int> const& vtx, vector<vector<int> > const& partsVert, vector<vector<int> > const& partsEdge, int const vTail, int const vHead, int& pProc, int& pIdx) {
  pProc = findVertIndex(vtx, vHead);
  int vLocIdx = vHead - vtx[pProc+1-1];
  int headEdgesStart = partsVert[pProc][vLocIdx];
  int headEdgesEnd = partsVert[pProc][vLocIdx+1];
  vector<int> headEdges(&partsEdge[pProc][headEdgesStart], &partsEdge[pProc][headEdgesEnd] + 1);
  const int* pIdxPtr = std::find(&partsEdge[pProc][headEdgesStart], &partsEdge[pProc][headEdgesEnd] + 1, vTail);
  pIdx = pIdxPtr - &partsEdge[pProc][0];
}

vector<vector<int> > genSymEdgeWgtsForParts(vector<int>& vtx, vector<vector<int> >& partsVert, vector<vector<int> >& partsEdge, int base=WGT_BASE) {
  // create constant weights
  vector<vector<int> > ret = genConWgtsForParts(partsEdge, -1);

  // go through all edges
  int iEdges = 0;
  for (vector<vector<int> >::iterator outerIt=ret.begin(); outerIt!=ret.end(); ++outerIt) {
    for (vector<int>::iterator innerIt=outerIt->begin(); innerIt!=outerIt->end(); ++innerIt) {
      // if not set yet
      if (*innerIt == -1) {
        // get edge coords
        int eProc = outerIt - ret.begin();
        int eIdx = innerIt - outerIt->begin();

        // get head and tail vertex indices
        int vTail = getTailVertexIndex(vtx, partsVert, partsEdge, eProc, eIdx);
        int vHead = partsEdge[eProc][eIdx];

        // current edge weight
        int curWeight = base + FAC_SYM * iEdges;

        // set edge weight
        ret[eProc][eIdx] = curWeight;

        // find pairing edge
        int pProc = -1; // TODO
        int pIdx = -1; // TODO
        findPairingEdgeIndices(vtx, partsVert, partsEdge, vTail, vHead, pProc, pIdx);

        // set pairing edge weight
        ret[pProc][pIdx] = curWeight;

        // increment edge count
        ++iEdges;
      }
    }
  }
  return ret;
}

void printGraphDataFile(vector<int>& vtx, vector<vector<int> >& partsVert, vector<vector<int> >& partsEdge) {
  string defType = "";
  defType.append(ARRAY_TYPE_DEF);
  defType.append(" ");
  defType.append(STRUCT_NAME);
  defType.append("::");

  cout << "// example graph data for " << vtx.size()-1 << " procs" << endl;
  // struct begin
  cout << "struct " << STRUCT_NAME << " {" << endl;

  // declare
  cout << ARRAY_TYPE_DECLARE << " "  << vecInfoSizeOnly<int>(vtx, "exVertglbtabRaw") << ";" << endl;
  cout << endl;

  cout << ARRAY_TYPE_DECLARE << " "  << vecInfoSizeOnly<int>(concatVecVec<int>(partsVert), "exVertloctabRaw") << ";" << endl;
  cout << ARRAY_TYPE_DECLARE << " "  << vecInfoSizeOnly<int>(concatVecVec<int>(genWgtsForParts(partsVert, WGT_BASE+MOD_WGTCMP)), "exVertwgtcmpRaw") << ";" << endl;
  cout << ARRAY_TYPE_DECLARE << " "  << vecInfoSizeOnly<int>(concatVecVec<int>(genWgtsForParts(partsVert, WGT_BASE+MOD_WGTVIS)), "exVertwgtvisRaw") << ";" << endl;
  cout << endl;

  // WARNING: genConWgtsForParts is used only because solely vector size matters here!
  cout << ARRAY_TYPE_DECLARE << " "  << vecInfoSizeOnly<int>(concatVecVec<int>(partsEdge), "exEdgeloctabRaw") << ";" << endl;
  cout << ARRAY_TYPE_DECLARE << " "  << vecInfoSizeOnly<int>(concatVecVec<int>(genConWgtsForParts(partsEdge, WGT_BASE+MOD_WGTCMP)), "exEdgewgtcmpRaw") << ";" << endl;
  cout << ARRAY_TYPE_DECLARE << " "  << vecInfoSizeOnly<int>(concatVecVec<int>(genConWgtsForParts(partsEdge, WGT_BASE+MOD_WGTVIS)), "exEdgewgtvisRaw") << ";" << endl;

  // struct end
  cout << "};" << endl << endl;

  // define
  cout << defType << vecInfo<int>(vtx, "exVertglbtabRaw") << ";" << endl;
  cout << endl;

  cout << defType << vecInfo<int>(concatVecVec<int>(partsVert), "exVertloctabRaw") << ";" << endl;
  cout << defType << vecInfo<int>(concatVecVec<int>(genWgtsForParts(partsVert, WGT_BASE+MOD_WGTCMP)), "exVertwgtcmpRaw") << ";" << endl;
  cout << defType << vecInfo<int>(concatVecVec<int>(genWgtsForParts(partsVert, WGT_BASE+MOD_WGTVIS)), "exVertwgtvisRaw") << ";" << endl;
  cout << endl;

  cout << defType << vecInfo<int>(concatVecVec<int>(partsEdge), "exEdgeloctabRaw") << ";" << endl;
  cout << defType << vecInfo<int>(concatVecVec<int>(genSymEdgeWgtsForParts(vtx, partsVert, partsEdge, WGT_BASE+MOD_WGTCMP)), "exEdgewgtcmpRaw") << ";" << endl;
  cout << defType << vecInfo<int>(concatVecVec<int>(genSymEdgeWgtsForParts(vtx, partsVert, partsEdge, WGT_BASE+MOD_WGTVIS)), "exEdgewgtvisRaw") << ";" << endl;
}

int main(int argc, char** argv) {

  int w = 5;
  int h = 3;

  int p = 3;

  bool verbose = false;
  
  // get possible arguments
  vector<string> args(argv + 1, argv + argc);
  vector<string>::iterator it;
  it=find(args.begin(), args.end(), "-v"); if (it!=args.end()) { verbose = true;}
  it=find(args.begin(), args.end(), "-?"); if (it!=args.end()) { cout << "Graph generator for a rectangular graph. Options: -w [width[=5]] -h [height[=3]] -p [procs[=3]] -v" << endl; exit(0);}
  it=find(args.begin(), args.end(), "--help"); if (it!=args.end()) { cout << "Graph generator for a rectangular graph. Options: -w [width[=5]] -h [height[=3]] -p [procs[=3]] -v" << endl; exit(0);}
  it=find(args.begin(), args.end(), "-w"); if (it!=args.end() && ++it!=args.end()) { w = strToNbr<int>(*it); }
  it=find(args.begin(), args.end(), "-h"); if (it!=args.end() && ++it!=args.end()) { h = strToNbr<int>(*it); }
  it=find(args.begin(), args.end(), "-p"); if (it!=args.end() && ++it!=args.end()) { p = strToNbr<int>(*it); }
 
  if (w<2 || h<2 || p<1) { cout << "Use reasonable width, height and parts values (current: w=" << w << ", h=" << h << ", p=" << p << ")." << endl; exit(-1); }
  
  vector<int> vertglb;
  vector<int> edgeglb;

  vector<vector<int> > partsVert;
  vector<vector<int> > partsEdge;

  genGlobal(w, h, vertglb, edgeglb);

  if (verbose) {
    cout << "serial:" << endl;
    cout << vecInfo<int>(vertglb, "verttab") << endl;
    cout << vecInfo<int>(edgeglb, "edgetab") << endl;
    cout << endl;
  }

  if (verbose) { cout << "parallel:" << endl; }
  vector<int> vtx(1, 0);
  for (int k=0; k<p; ++k) {
    vector<int> vertloc;
    vector<int> edgeloc;

    getPart(p, k, vertglb, edgeglb, vertloc, edgeloc);

    vtx.push_back(vertloc.size()-1+vtx[k]);

    if (verbose) {
      cout << "[" << k << "]:" << endl;
      cout << vecInfo<int>(vertloc, "vertloctab") << endl;
      cout << vecInfo<int>(edgeloc, "edgeloctab") << endl;
      cout << endl;
    }

    partsVert.push_back(vertloc);
    partsEdge.push_back(edgeloc);
  }
  if (verbose) { cout << vecInfo<int>(vtx, "vtx") << endl; }

  printGraphDataFile(vtx, partsVert, partsEdge);
}
