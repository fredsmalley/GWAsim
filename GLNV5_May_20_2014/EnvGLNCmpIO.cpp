// EnvGLNCmpIO.cpp
//
// Joe Song
// Created: December 8, 2009.  
//
// Updated: 
//
//   April 7, 2011. Joe Song. Updated saveIntxCmp() to the latest version
//     developed on March 14, 2010.
//
//   December 4, 2011. MS. Updated EnvGLNCmp::saveIntxCmp().
//     - Changed name from EnvGLNCmp::saveIntxCmp() to 
//       EnvGLNCmp::save(const string & file)
//     - Changed header from "SameGTT\tSameParents\tTransProbTableDistance"
//       to "Type\tSameParents\tp_d"
//     - Include p_t and p_z as the last two columns in the output file.
//   
//   Jan 11, 2012. MS. Automatically determine to add double quote when
//     generating DOT file.  Redundant double quotes cause problem in DOT
//     to PDF file conversion.
//
//   Feb 26, 2012. MS. Added void EnvGLNCmp::saveTopology(const string & 
//     fileTopology) to save all non-null interactions where both child and 
//     parents have changed working zones.  There is another version to save 
//     selected type of comparative interactions.
//
//   June 22, 2013. MS. Now the EnvGLNCmp::saveTopology() output all parents
//     including all conditions for each child with duplicates removed.
//
//   May 3, 2014. MS. Amended quote_as_needed() to remove all intermediate
//     double quote symbols within the input string.
//
//   May 6, 2014. MS. Added member functions saveDOTFilebyCondition() and
//     saveDOTCondition() to support more than two conditions for generating
//     the output graph file in DOT format. This function thus replaces
//     saveDotFilek2() which supports exactly two conditions.

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <algorithm>
using std::find;
using std::sort;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;

#include "SetOps.h"
#include "GLNGlobals.h"
#include "EnvGLNCmp.h"

void EnvGLNCmp::save(const string & file) const 
// Result of comparative analysis: 
//   1. used for performance evaluation
//   2. used also for further analysis 
{
    if(! file.empty() ) {
        
		ofstream ofs(file.c_str());
		
        // Header
        
		// ofs << "SameGTT\tSameParents\tTransProbTableDistance"; MS Dec 4, 2011
        
        ofs << "Type\tSameParents\tp_d";  // MS Dec 4, 2011
        
        for (size_t k=0; k<m_TransTables[0].size(); k++) {
            ofs << "\t";
            ofs << "Parents" << k+1;
        }
        
        ofs << "\tCommonParents\tp_c\tp_t\tp_z(child)\tp_z(parents)\tChildName" << endl;
        
        
        // Body
        
		for(size_t i=0; i<m_diffTransTables.size(); i++) {
			
			if(true) { // Two comparative interaction types
				
				// if(m_diffTransTables[i].getpValue() <= m_gc.m_alpha) {
				if(m_types[i] == REL_DIFF || m_types[i] == ABS_DIFF) {
					
					ofs << "D";
					
					// } else if(m_pooledTransTables[i].getpValue() <= m_gc.m_alpha) {
				} else if(m_types[i] == CONSERVED || m_types[i] == NULL_INTX) {
					
					ofs << "C"; 
					
				}
				
			} else { // Four comparative interaction types
				
				// if(m_diffTransTables[i].getpValue() <= m_gc.m_alpha) {
				if(m_types[i] == REL_DIFF || m_types[i] == ABS_DIFF) {
					
					if(m_pooledTransTables[i].getpValue() <= m_gc.m_alpha) {
						ofs << "B"; // The interaction has both significant 
						// heterogeneous and homoegeneous components
					} else {
						ofs << "D"; // The interaction has only significant 
						// heterogeneous components 
					}
					
					// } else if(m_pooledTransTables[i].getpValue() <= m_gc.m_alpha) {
				} else if(m_types[i] == CONSERVED) {
					
					ofs << "C"; // The interaction has a significant 
					// homoegeneous components.
					
				} else if (m_types[i] == NULL_INTX) {
					
					ofs << "N"; // Null interaction.
					
				}
			}
			
			ofs << '\t';
			ofs << '-' << '\t'; // SameParents entry is not used.  Reserved for 
                // reconstruct-then-compare
            
			ofs << m_diffTransTables[i].getpValue();
            
            for(size_t k=0; k<m_TransTables[i].size(); k++) {
                
                ofs << '\t';                        
                
                for(size_t j=0; j < m_TransTables[i][k].getParents().size(); 
                    j++) { // Joe Song. March 12, 2010
					
                    if(j > 0) ofs << ',';
					
                    ofs << m_TransTables[i][k].getParents()[j];  // Joe Song. March 12, 2010
                }
            }
			
			ofs << '\t';
            
            for(size_t j=0; j < m_pooledTransTables[i].getParents().size(); j++) {
                
                if(j > 0) ofs << ',';
                
                ofs << m_pooledTransTables[i].getParents()[j]; 
            }
            
            ofs << '\t' << m_pooledTransTables[i].getpValue();
			
            ofs << '\t' << m_pchisqs_t[i];
            
            ofs << '\t' << m_childpchisqs_z[i];

            ofs << '\t' << m_pchisqs_z[i];
            
            ofs << '\t' << m_glnPooled.getNodeName(i);
            
			ofs << endl;
		}
		
		ofs.close();
	}
}

void EnvGLNCmp::print(const GeneralizedLogicalNetwork & gln, 
					  const vector<TransitionTable> & tts) const
{
	cout << endl << endl;

	cout << "id\t";
	cout << "name\t";
	cout << "type\t";
	cout << "num.parents\t";
	cout << "parents\t";
	cout << "offsets\t";
	cout << "p.value\t";
	cout << "chisq\t";
	cout << "df" << endl;

	for(size_t i=0; i < gln.size(); i++){

		cout << i+1 << "\t";
		cout << gln.getNodeName(i) << "\t";
		cout << gln.getNodeType(i) << "\t";

		tts[i].heading();

	}

	// cout << "Reconstructed generalized logical network p-value = " << m_p_val << endl;

	print_adjusted(gln, tts, m_gc.m_alpha);	

	for(size_t i=0; i < gln.size(); i++){

		cout << "Significant transition table for node " << i+1 <<": ";

		if(tts[i].nrow()==0) {
			cout << "None." << endl;
		} else {
			cout << endl;
			tts[i].print();
		}
	}

	gln.printGTTs();

}

void EnvGLNCmp::print() const
{
	cout << endl << endl;
	cout << "*************************************************************" << endl;
	cout << "Significant comparative interactions identified:" << endl;
	cout << "___________________________________________________________________________" << endl;
    /*
	cout << "id\t";
	cout << "name\t";
	//cout << "node.type\t";
	//cout << "cmp.intx.type\t";
	//cout << "num.parents\t";
	cout << "parents\t";
	//cout << "offsets\t";
	cout << "chisq_d\t";
	cout << "df_d\t";
	cout << "p_d\t";
	cout << "chisq_c\t";
	cout << "df_c\t";
	cout << "p_c\t";
	cout << "chisq_t\t";
	cout << "df_t\t";
	cout << "p_t\t";
	cout << "chisq_z\t";
	cout << "df_z\t";
	cout << "p_z\t";
	cout << endl;
    */
    
	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == ABS_DIFF && m_pooledTransTables[i].getpValue() > m_gc.m_alpha
			&& m_childWorkingZones[i] == CHANGED) {
            
            if (total == 0) {
                cout << "_________________________________________________________________" << endl;
                cout << "ABSOLUTE DIFFERENTIAL, without significant homogeneous component:" << endl;
            }            
			print(i);
            total ++;
		}
	}
	
	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == ABS_DIFF && m_pooledTransTables[i].getpValue() <= m_gc.m_alpha
			&& m_childWorkingZones[i] == CHANGED) {

            if (total == 0) {
                cout << endl;
                cout << "______________________________________________________________" << endl;
                cout << "ABSOLUTE DIFFERENTIAL, with significant homogeneous component:" << endl;
            }            
			print(i);
            total ++;

		} else {
            // cout << "DEBUG (SHALL NOT PRINT) ...";
            // print(i);
        }
        
	}

	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == REL_DIFF && m_pooledTransTables[i].getpValue() > m_gc.m_alpha
			&& m_childWorkingZones[i] == CHANGED) {

            if (total == 0) {
                cout << endl;
                cout << "_________________________________________________________________" << endl;
                cout << "RELATIVE DIFFERENTIAL, without significant homogeneous component:" << endl;
            }            
			print(i);
            total ++;
		}
	}
	
	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == REL_DIFF && m_pooledTransTables[i].getpValue() <= m_gc.m_alpha
			&& m_childWorkingZones[i] == CHANGED) {

            if (total == 0) {
                cout << endl;
                cout << "______________________________________________________________" << endl;
                cout << "RELATIVE DIFFERENTIAL, with significant homogeneous component:" << endl;
            }            
			print(i);
            total ++;
		}
	}

	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == CONSERVED && m_workingZones[i] == CHANGED 
           && m_childWorkingZones[i] == CHANGED) {

            if (total == 0) {
                cout << endl;
                cout << "________________________________________" << endl;
                cout << "CONSERVED, but with WORKING ZONE CHANGE:" << endl;
            }            
			print(i);
            total ++;
		}
	}

	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == CONSERVED && m_workingZones[i] == UNCHANGED 
			&& m_pchisqs_z[i] > 1 - m_gc.m_alpha
			&& m_childWorkingZones[i] == UNCHANGED  
			&& m_childpchisqs_z[i] > 1 - m_gc.m_alpha) {

            if (total == 0) {
                cout << endl;
                cout << "___________________________________________________" << endl;
                cout << "CONSERVED, without working zone change: (IDENTICAL)" << endl;
            }            
			print(i);
            total ++;
		}
	}
	cout << endl;

	cout << endl;

}

//currently screen print is different for k=1 and k=2
void EnvGLNCmp::printk2() const
{
	cout << endl << endl;
	cout << "*************************************************************" << endl;
	cout << "Significant comparative interactions identified:" << endl;
	cout << "___________________________________________________________________________" << endl;
    
	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == ABS_DIFF && m_pooledTransTables[i].getpValue() > m_gc.m_alpha) {
            
            if (total == 0) {
                cout << "_________________________________________________________________" << endl;
                cout << "ABSOLUTE DIFFERENTIAL, without significant homogeneous component:" << endl;
            }            
			print(i);
            total ++;
		}
	}
	
	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == ABS_DIFF && m_pooledTransTables[i].getpValue() <= m_gc.m_alpha) {

            if (total == 0) {
                cout << endl;
                cout << "______________________________________________________________" << endl;
                cout << "ABSOLUTE DIFFERENTIAL, with significant homogeneous component:" << endl;
            }            
			print(i);
            total ++;

		}
	}

	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == REL_DIFF && m_pooledTransTables[i].getpValue() > m_gc.m_alpha) {

            if (total == 0) {
                cout << endl;
                cout << "_________________________________________________________________" << endl;
                cout << "RELATIVE DIFFERENTIAL, without significant homogeneous component:" << endl;
            }            
			print(i);
            total ++;
		}
	}
	
	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == REL_DIFF && m_pooledTransTables[i].getpValue() <= m_gc.m_alpha) {

            if (total == 0) {
                cout << endl;
                cout << "______________________________________________________________" << endl;
                cout << "RELATIVE DIFFERENTIAL, with significant homogeneous component:" << endl;
            }            
			print(i);
            total ++;
		}
	}

	for(int i=0, total=0; i < (int) m_glnPooled.size(); i++){
		if(m_types[i] == CONSERVED) {

            if (total == 0) {
                cout << endl;
                cout << "___________" << endl;
                cout << "CONSERVED:" << endl;
            }            
			print(i);
            total ++;
		}
	}

	cout << endl;

	cout << endl;

}

void EnvGLNCmp::print(int i) const
{
    cout << i+1 << "\t";
    cout << m_glnPooled.getNodeName(i) << "\t";
    
    // cout << m_glnPooled.getNodeType(i) << "\t";
    //cout << m_types[i] << "\t";
    
    /*
    vector<int>parentsPooled = m_pooledTransTables[i].getParents();
    // vector<int>parentsDiff = m_diffTransTables[i].getParents();
    
    for(size_t k=0; k<parentsPooled.size(); k++) {
        cout <<	m_glnPooled.getNodeName(parentsPooled[k]-1) 
        // << m_glnDiff.getNodes()[parentsDiff[k]-1].getName() 
        << ",";
    }
     */

    cout << endl;
    
    cout << showpoint;
    
    // cout.width(8);
    
    cout << "chi2d=" << m_diffTransTables[i].getChisq() << "\t" 
    << "vd=" << m_diffTransTables[i].getDf() << "\t" 
    << "pd=" << m_diffTransTables[i].getpValue() << "\t(heterogeneity)" << endl;
    
    cout << "chi2c=" << m_pooledTransTables[i].getChisq() << "\t" 
    << "vc=" << m_pooledTransTables[i].getDf() << "\t" 
    << "pc=" << m_pooledTransTables[i].getpValue() << "\t(homogeneity)" << endl;
    
    cout << "chi2t=" << m_chisqs_t[i] << "\t" 
    << "vt=" << m_dfs_t[i] << "\t" 
    << "pt=" << m_pchisqs_t[i] << "\t(total strength)" << endl;

    cout << "chi2z=" << m_childchisqs_z[i] << "\t" 
    << "vz=" << m_childdfs_z[i] << "\t" 
    << "pz=" << m_childpchisqs_z[i] << "\t(child working zone change)" << endl;
    
    cout << "chi2z=" << m_chisqs_z[i] << "\t" 
    << "vz=" << m_dfs_z[i] << "\t" 
    << "pz=" << m_pchisqs_z[i] << "\t(parent working zone change)" << endl;
        
    // print contingency tables

    for(size_t j=0; j<m_TransTables[i].size(); j++) {
        cout << "... Contingency table under condition " << j+1 << ":" << endl;
        m_TransTables[i][j].print();
    }
    
    cout << "... Pooled contingency table of all conditions:" << endl;

    m_pooledTransTables[i].print();
    
    cout << "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~" << endl;
}


/*
void EnvGLNCmp::print() const
{
	cout<<endl;
	if (checkdiff(m_diffTransTables,m_pooledTransTables,m_gc.m_alpha)==2)
	{
		//FILE * outputTrj;
		//outputTrj = fopen("compresult.txt","a+");
		//fprintf(outputTrj,"%d\n", 1);
		//fclose(outputTrj);
		cout<<"The input trajectories imply different networks."<<endl;
		//cout<<"The input trajectories imply different network."<<endl;
	}
	else
	{
		cout<<"The input trajectories imply the same network."<<endl;
	}
	print(m_glnPooled, m_pooledTransTables);
	print(m_glnDiff, m_diffTransTables);
}
*/

string quote_as_needed(const string & str)
{ // apply double quotes to str only if str is not already in quote.
    
    string processed;
            
    size_t pos = find(str.begin(), str.end(), '\"') - str.begin();
    
    if( pos == str.size() ) {
        
        processed = "\"" + str + "\"";
        
    } else if( str[0] == '\"' && * str.rbegin() == '\"' ) {
        
        processed = str;
        
    } else if( str[0] == '\"' && * str.rbegin() != '\"' ) {
       
        processed = str;
        
        size_t i = 1;
        while (i < str.size()) {
            size_t j = find(str.begin()+i, str.end(), '\"') - str.begin();
            if (j < str.size()) {
                processed[j] = '.'; // replace double quote by a dot
            }
            i = j+1;
        }

        processed += "\"";

    } else {
        cerr << "ERROR: Invalid string for dot in " << str << "!" << endl;
        GLNExit(EXIT_FAILURE);
    }
    
    return processed;
}

void EnvGLNCmp::saveDotInteraction(const string & filename)
{
	ofstream out(filename.c_str()); 
	//green,blue,chocolate represents common,absolute and relative differential respectively
	// string colors[3] = {"green","blue","chocolate"};
	string colors[3] = {"palegreen","steelblue1","aliceblue"};
    
	out.precision(2);
    
	out << "digraph GRN {" << endl
    << "ratio=auto;" << endl
    << "margin=0;" <<endl
    << "node[fontname=\"Helvetica\"];" << endl
    << "node[fontcolor=\"black\"];" << endl
    << "rankdir=LR;" << endl
    << "labelloc=t;" << endl;
    
	for(int i=0; i < (int) getGLNSize(); i++)
	{
		// int dest = (int) i+1;
		//0 denotes null interation, 1 denotes relative differential, 2 denotes absolute differential, 3 denotes homogenous interaction
		if(getType(i)==REL_DIFF && m_childWorkingZones[i] == CHANGED) // 1)	//relative differential
		{
			vector<int> delays = getDiffGLN().getGTT(i).getDelays();
			string nodeNameDest = quote_as_needed(getDiffGLN().getNodeName(i));
			for(size_t j=0; j<getDiffTransTable(i).getParents().size(); j++)
			{
                int src = getDiffTransTable(i).getParents()[j];
                string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1));
                
                // out	<< "\"" << nodeNameSrc << "\" " << "-> " << "\"" << nodeNameDest << "\"";
                out	<< nodeNameSrc << " -> " << nodeNameDest;
                
                //out << " [label=" << delays[j] <<","  
                out << " [label=\"" << delays[j] <<";"<< getDiffTransTable(i).getpValue() <<"\","   
				<< "arrowhead=normal,penwidth=3,color="<<colors[2]<<",style=solid, ];" << endl;
			}
		}
		else if(getType(i)==ABS_DIFF && m_childWorkingZones[i] == CHANGED) // 2)   //absolute differential
		{
			vector<int> delays = getDiffGLN().getGTT(i).getDelays();
			string nodeNameDest = quote_as_needed(getDiffGLN().getNodeName(i));
			for(size_t j=0; j<getDiffTransTable(i).getParents().size(); j++)
			{
                int src = getDiffTransTable(i).getParents()[j];
                string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1));
                out	<< nodeNameSrc << " -> " << nodeNameDest;
                //out << " [label=" << delays[j] <<","  
                out << " [label=\"" << delays[j] <<";"<< getDiffTransTable(i).getpValue() <<"\","  
				<< "arrowhead=normal,penwidth=3,color="<<colors[1]<<",style=solid, ];" << endl;
			}
		}
		else if(getType(i)==CONSERVED && m_workingZones[i] == CHANGED && m_childWorkingZones[i] == CHANGED) // 3)   //homogenous interaction
		{
			vector<int> delays = getPooledGLN().getGTT(i).getDelays();
			string nodeNameDest = quote_as_needed(getPooledGLN().getNodeName(i));
			for(size_t j=0; j<getPooledTransTable(i).getParents().size(); j++)
			{
                int src = getPooledTransTable(i).getParents()[j];
                string nodeNameSrc = quote_as_needed(getPooledGLN().getNodeName(src-1));
                out	<< nodeNameSrc << " -> " << nodeNameDest;
                //out << " [label=" << delays[j] <<","  
                out << " [label=\"" << delays[j] <<";"<< getPooledTransTable(i).getpValue() <<"\","  
				<< "arrowhead=normal,color="<<colors[0]<<",style=solid, ];" << endl;
			}
		}
	}
	out << endl << endl;
    
	for(int i=0; i < (int) getGLNSize(); i++){
        
        // int dest = (int) i+1;
        string nodeName = quote_as_needed(getPooledGLN().getNodeName(i));
        
        //if(getType(i) != NULL_INTX ) { // > 0 ) {
        //}
        
        switch(getType(i)) {
			case NULL_INTX: // 0:
				// out << "[shape=ellipse,style=filled,fillcolor=white];";
				break;
			case REL_DIFF: // 1:
				if(m_childWorkingZones[i] == CHANGED) {
					out << nodeName;
					out << "[shape=ellipse,style=filled,fillcolor="<<colors[2]<<"];";
					out << endl;
				}
				break;
			case ABS_DIFF: // 2:
				if(m_childWorkingZones[i] == CHANGED) {
					out << nodeName;
					out << "[shape=ellipse,style=filled,fillcolor="<<colors[1]<<"];";
					out << endl;
				}
				break;
			case CONSERVED: // 3:
				if(m_workingZones[i] == CHANGED && m_childWorkingZones[i] == CHANGED) {
					out << nodeName;
					out<<  "[shape=ellipse,style=filled,fillcolor="<<colors[0]<<"];";
					out << endl;
                    
				}
        }
	}
	out << endl << "}";
	out.close();
}

int choose_penwidth(double pval)
{
    int pw=0;

    if(pval <= 0.001) 
        pw = 7;
    else if( 0.001 < pval && pval <= 0.01 ) 
        pw = 5;
    else if( 0.01 < pval && pval <= 0.05 ) 
        pw = 5;
    else if( 0.05 < pval && pval <= 0.1 )
        pw = 4;
    else if( 0.1 < pval && pval <= 0.2 )
        pw = 3;
    else if( 0.2 < pval && pval <= 0.5 )
        pw = 2;
    else if( 0.5 < pval )
        pw = 1;

    return pw;
}

// does not support differential parents, modified in saveDotFileK2, changed by Yang Zhang 10.24.2012
void EnvGLNCmp::saveDotFile(const string & filename)
{
	vector <int> subgraphAbsDiff;
	vector <int> subgraphRelDiff;
	vector <int> subgraphConservedZoneChanged;
	vector <int> subgraphIdentical;
	vector <int> subgraphNullIntx;

	for(int i=0; i < (int) getGLNSize(); i++) {
		
		if(getType(i)==REL_DIFF && m_childWorkingZones[i] == CHANGED) { //relative differential

			subgraphRelDiff.push_back(i);
			
		} else if(getType(i)==ABS_DIFF && m_childWorkingZones[i] == CHANGED) { //absolute differential
		
			subgraphAbsDiff.push_back(i);
			
		} else if(getType(i)==CONSERVED 
				  && m_workingZones[i] == CHANGED 
				  && m_childWorkingZones[i] == CHANGED) { // conserved with working zone change

			subgraphConservedZoneChanged.push_back(i);

		} else if(getType(i)==CONSERVED 
				  && m_workingZones[i] == UNCHANGED 
				  && m_pchisqs_z[i] > 1 - m_gc.m_alpha
				  && m_childWorkingZones[i] == UNCHANGED
		     	  && m_childpchisqs_z[i] > 1 - m_gc.m_alpha ) { // conserved with working zone change

			subgraphIdentical.push_back(i);

		} else {

			subgraphNullIntx.push_back(i);

		}
	}

	ofstream out(filename.c_str());
	out.precision(2);

	//green,blue,chocolate represents common,absolute and relative differential respectively
	// string colors[3] = {"green","blue","chocolate"};

	string colors[] = {"palegreen", "steelblue1", "aliceblue", "forestgreen", "gray"};

	time_t rawtime;
	time ( &rawtime );
	string strTime = ctime(&rawtime);
	strTime.erase(strTime.end()-1);

	out << "digraph GLNComparativeModeling {" << endl
		<< "label=\"GLN Comparative Modeling: " << filename << "\t" << strTime << "\"" << endl
		<< "ratio=auto;" << endl
		<< "margin=0;" <<endl
        << "fontname=\"Helvetica\";" << endl
        << "edge[fontname=\"Helvetica\"];" << endl
		<< "node[fontname=\"Helvetica\"];" << endl
        << "node[fontcolor=\"black\"];" << endl
		<< "rankdir=LR;" << endl
		<< "labelloc=t;" << endl;

	out << "subgraph \"cluster_REL_DIFF_WORKINGZONE_CHANGED\" {" << endl;
	out << "label=\"" << "Relatively Differential (Working Zone Changed)\";" << endl;

	for(int n=0; n < (int) subgraphRelDiff.size(); n++)
	{
		int i = subgraphRelDiff[n];
		// int dest = (int) i+1;

		// vector<int> delays = getDiffGLN().getGTT(i).getDelays();
        const vector<int> & delays = getTransTables(i)[0].getDelays();

		string nodeNameDest = quote_as_needed(getDiffGLN().getNodeName(i));

		// for(size_t j=0; j<getDiffTransTable(i).getParents().size(); j++) {
		for(size_t j=0; j<getTransTables(i)[0].getParents().size(); j++) {
            
            // int src = getDiffTransTable(i).getParents()[j];
            int src = getTransTables(i)[0].getParents()[j];
            
			string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1));
			out	<< nodeNameSrc << " -> " << nodeNameDest;
			//out << " [label=" << delays[j] <<","  
            
            int penwidth = choose_penwidth( getDiffTransTable(i).getpValue() );
            
			out << " [label=\"" << delays[j] <<";"<< getDiffTransTable(i).getpValue() <<"\","   
				<< "arrowhead=normal,penwidth=" << penwidth 
                << ",color="<<colors[2]<<",style=solid, ];" << endl;
		}

		out << nodeNameDest;
		out << "[shape=ellipse,style=filled,fillcolor="<<colors[2]<<"];";
		out << endl;

	}

	out << "}" << endl << endl;

	out << "subgraph \"cluster_ABS_DIFF_WORKINGZONE_CHANGED\" {" << endl;
	out << "label=\"" << "Absolutely Differential (Working Zone Changed)\";" << endl;

	for(int n=0; n < (int) subgraphAbsDiff.size(); n++) {

		int i = subgraphAbsDiff[n];
		// int dest = (int) i+1;

		// vector<int> delays = getDiffGLN().getGTT(i).getDelays();
        const vector<int> & delays = getTransTables(i)[0].getDelays();
        
        string nodeNameDest = quote_as_needed(getDiffGLN().getNodeName(i));
        
		// for(size_t j=0; j<getDiffTransTable(i).getParents().size(); j++)
        for(size_t j=0; j<getTransTables(i)[0].getParents().size(); j++)
		{
			// int src = getDiffTransTable(i).getParents()[j];
            int src = getTransTables(i)[0].getParents()[j];
            
			string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1));
			out	<< nodeNameSrc << " -> " << nodeNameDest;
			//out << " [label=" << delays[j] <<","  

            int penwidth = choose_penwidth( getDiffTransTable(i).getpValue() );

			out << " [label=\"" << delays[j] <<";"<< getDiffTransTable(i).getpValue() <<"\","  
            << "arrowhead=normal,penwidth=" << penwidth 
            << ",color="<<colors[1]<<",style=solid, ];" << endl;
		}
		out << nodeNameDest;
		out << "[shape=ellipse,fontcolor=white,style=filled,fillcolor="<<colors[1]<<"];";
		out << endl;

	}
	out << "}" << endl << endl;

	out << "subgraph \"cluster_CONSERVED_WORKINGZONE_CHANGED\" {" << endl;
	out << "label=\"" << "Conserved (Working Zone Changed)\";" << endl;

	for(int n=0; n < (int) subgraphConservedZoneChanged.size(); n++) {

		int i = subgraphConservedZoneChanged[n];
		// int dest = (int) i+1;

		vector<int> delays = getPooledGLN().getGTT(i).getDelays();
		string nodeNameDest = quote_as_needed(getPooledGLN().getNodeName(i));
		for(size_t j=0; j<getPooledTransTable(i).getParents().size(); j++)
		{
			int src = getPooledTransTable(i).getParents()[j];
			string nodeNameSrc = quote_as_needed(getPooledGLN().getNodeName(src-1));
			out	<< nodeNameSrc << " -> " << nodeNameDest;
			//out << " [label=" << delays[j] <<","  
            
            int penwidth = choose_penwidth( getPooledTransTable(i).getpValue() );

			out << " [label=\"" << delays[j] <<";"<< getPooledTransTable(i).getpValue() <<"\","  
            << "arrowhead=normal,penwidth=" << penwidth 
            << ",color="<<colors[0]<<",style=solid, ];" << endl;
		}

		out << nodeNameDest;
		out<<  "[shape=ellipse,style=filled,fillcolor="<<colors[0]<<"];";
		out << endl;

	}

	out << "}" << endl << endl;

	out << "subgraph \"cluster_IDENTICAL\" {" << endl;
	out << "label=\"" << "Conserved (Working Zone Unchanged)\";" << endl;

	for(int n=0; n < (int) subgraphIdentical.size(); n++) {

		int i = subgraphIdentical[n];
		// int dest = (int) i+1;

		vector<int> delays = getPooledGLN().getGTT(i).getDelays();
		string nodeNameDest = quote_as_needed(getPooledGLN().getNodeName(i));
		for(size_t j=0; j<getPooledTransTable(i).getParents().size(); j++)
		{
			int src = getPooledTransTable(i).getParents()[j];
			string nodeNameSrc = quote_as_needed(getPooledGLN().getNodeName(src-1));
			out	<< nodeNameSrc << " -> " << nodeNameDest;

            int penwidth = choose_penwidth(getPooledTransTable(i).getpValue());

			out << " [label=\"" << delays[j] <<";"<< getPooledTransTable(i).getpValue() <<"\","  
            << "arrowhead=normal,penwidth=" << penwidth 
            << ",color="<<colors[3]<<",style=solid];" << endl;
		}

		out << nodeNameDest;
		out<<  "[shape=ellipse,style=filled,fontcolor=white,fillcolor="<<colors[3]<<"];";
		out << endl;

	}
	out << "}" << endl << endl;

	// out << "subgraph \"cluster_IDENTICAL\" {" << endl;
	// out << "label=\"" << "Conserved (Working Zone Unchanged)\";" << endl;
    
    if (false) { // show detected interactions that are insignificant relative to the threshold 
        for(int n=0; n < (int) subgraphNullIntx.size(); n++) {
            
            int i = subgraphNullIntx[n];
            
            if(getDiffTransTable(i).getParents().size() > 0) {
                
                // vector<int> delays = getDiffGLN().getGTT(i).getDelays();
                // const vector<int> & delays = getTransTables(i)[0].getDelays();

                string nodeNameDest = quote_as_needed(getDiffGLN().getNodeName(i));

                // for(size_t j=0; j < getDiffTransTable(i).getParents().size(); j++) {
                for(size_t j=0; j< getTransTables(i)[0].getParents().size(); j++) {
                    // int src = getDiffTransTable(i).getParents()[j];
                    int src = getTransTables(i)[0].getParents()[j];
                    string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1));
                    out	<< nodeNameSrc << " -> " << nodeNameDest;
                    
                    int penwidth = 1; // choose_penwidth( getDiffTransTable(i).getpValue() );
                    
                    out << " [" // label=\"" << delays[j] << ";" << getDiffTransTable(i).getpValue() << "\","   
                    << "arrowhead=normal,penwidth=" << penwidth 
                    << ",color=" << colors[4] << ",style=dashed];" << endl;
                }
                
                out << nodeNameDest;
                out << "[shape=ellipse,style=filled,fillcolor=" << colors[4] << "];";
                out << endl;
            } 
        }
        // out << "}" << endl << endl;
    }    
    
	out << endl << "}";
	out.close();
}

void EnvGLNCmp::saveDOTCondition(ofstream & out, size_t cond, const string colors[],
                                 bool showAllNodes, bool mergeNodes,
                                 vector<bool> & used, vector<bool> & isParent) const
{
    if(! mergeNodes) {
        out << "subgraph \"cluster_condition" << cond+1 << "\" {" << endl;
        out << "label=\"" << "condition " << cond+1 << "\";" << endl;
    }
    
    string appendix = mergeNodes ? "" : "C"+std::to_string(cond+1);
    
    const GeneralizedLogicalNetwork & diffGLN = getDiffGLN();
    
	for(size_t i=0; i < getGLNSize(); i++) {
        string nodeNameDest = quote_as_needed(diffGLN.getNodeName(i) + appendix);
        
        if((getType(i)==REL_DIFF||getType(i)==ABS_DIFF)) {
            //first condition
            for(size_t j=0; j<getTransTables(i)[cond].getParents().size(); j++) {
                int src = getTransTables(i)[cond].getParents()[j];
                string nodeNameSrc = quote_as_needed(
                                                     diffGLN.getNodeName(src-1)
                                                     +appendix
                                                     );
                out	<< nodeNameSrc << " -> " << nodeNameDest;
                
                int penwidth = choose_penwidth( getDiffTransTable(i).getpValue() );
                
                out << " [arrowhead=normal,penwidth=" << penwidth
                << ",color="<<colors[cond+1]<<",style=solid, ];" << endl;
                
                used[src-1] = true;  // the node is used for an edge
                isParent[src-1] = true;
            }
            
            used[i] = true;
            
        } else if(getType(i)==CONSERVED) {
            
            const vector<int> & parents = getTransTables(i)[cond].getParents();
            
            for(size_t j=0; j<parents.size(); j++) {
                int src = parents[j];
                string nodeNameSrc = quote_as_needed(diffGLN.getNodeName(src-1)
                                                     +appendix);
                out	<< nodeNameSrc << " -> " << nodeNameDest;
                
                int penwidth = choose_penwidth( getPooledTransTable(i).getpValue() );
                
                out << " [arrowhead=normal,penwidth=" << penwidth
                    << ",color="<<colors[0]<<",style=solid, ];" << endl;
                
                used[src-1] = true;  // the node is used for an edge
                isParent[src-1] = true;
            }
            
            used[i] = true;
            
        } /*else {  //as discussed with Joe, do not need to display null edges
           //first condition
           for(size_t j=0; j<getTransTables(i)[0].getParents().size(); j++)
           {
           int src = getTransTables(i)[0].getParents()[j];
           string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1)+"C1");
           out	<< nodeNameSrc << " -> " << nodeNameDest;
           
           int penwidth = 1;
           
           out << " [arrowhead=normal,penwidth=" << penwidth
           << ",color="<<colors[3]<<",style=solid, ];" << endl;
           }
           
           }*/
    }
    
    string colorScheme;
    if (! mergeNodes) {
        colorScheme = "fontcolor=white,color=transparent,style=filled,fillcolor="
            + colors[cond+1];
    } else {
        colorScheme = "fontcolor=black,color=" + colors[0];
    }
    
    if (!mergeNodes) {
        
        for(size_t i=0; i < getGLNSize(); i++) {
            
            if(used[i]) {
                
                const string & nodeName = diffGLN.getNodeName(i);
                string nodeLabel = quote_as_needed(nodeName+appendix);
                string shape = isParent[i] ? "diamond" : "ellipse";
                
                out << nodeLabel;
                out << "["
                << "label=" << quote_as_needed(nodeName)
                << ','
                << "shape=" << shape
                << ','
                << colorScheme
                << "];";
                out << endl;
                
                if(isParent[i]) {
                    out << "root ->" << nodeLabel << "[color=transparent];" << endl;
                }
            }
        }

        out << "root[style=invisible];" << endl;

    }

    
    if(! mergeNodes) {
        out << "}" << endl << endl;
    } else {
        out << endl;
    }
}

void EnvGLNCmp::saveDOTNodes(ofstream & out, const vector<bool> & used,
                                const vector<bool> & isParent) const
{
    const GeneralizedLogicalNetwork & diffGLN = getDiffGLN();
    
    string colorScheme = "fontcolor=black,color=forestgreen";

    for(size_t i=0; i < getGLNSize(); i++) {
        
        if(isParent[i]) {
            const string & nodeName = diffGLN.getNodeName(i);
            string nodeLabel = quote_as_needed(nodeName);
            out << "root ->" << nodeLabel << "[color=transparent];" << endl;
        }
    }
    
    for(size_t i=0; i < getGLNSize(); i++) {
        
        if(used[i]) {
            
            const string & nodeName = diffGLN.getNodeName(i);
            string nodeLabel = quote_as_needed(nodeName);
            string shape = isParent[i] ? "diamond" : "ellipse";
            
            out << nodeLabel;
            out << "["
            << "label=" << quote_as_needed(nodeName)
            << ','
            << "shape=" << shape
            << ','
            << colorScheme
            << "];";
            out << endl;
        }
    }
    
    out << "root[style=invisible];" << endl;

}

void EnvGLNCmp::saveDOTFilebyCondition(const string & filename) const
{ //support visualization with different parents

    bool showAllNodes = false; // Display all singleton nodes (true) or not (false)
    bool mergeNodes = true; // true: Use one unique node for the same variable across
                            // all conditions; false: each node is used by only
                            // one condition
    
	// colors[0] represents conserved interaction.
    // colors[1:11] represents condition 1 to 11.
	string colors[] = {"darkolivegreen4", "tomato", "royalblue", "forestgreen", "steelblue1", "chocolate", "gray", "salmon", "navy", "brown", "yellow", "red", "blue", "black", "palegreen"};
    
    if (getTransTables(0).size() > 14) {
        cerr << "WARNING: no DOT file is generated because only at most 14 conditions are supported in DOT graphs."
        << endl;
        return;
    }

    vector<bool> used(getGLNSize(), showAllNodes);
    vector<bool> isParent(getGLNSize(), false);

	ofstream out(filename.c_str());
	out.precision(2);
    
	out << "digraph GLNComparativeModeling {" << endl
    << "ratio=auto;" << endl
    << "margin=0;" <<endl
    << "fontname=\"Helvetica\";" << endl
    << "edge[fontname=\"Helvetica\",len=2];" << endl
    << "node[fontname=\"Helvetica\"];" << endl
    << "node[fontcolor=\"black\"];" << endl
    << "rankdir=TB;" << endl
    << "outputorder=edgesfirst;" << endl
    // << "overlap=prism;" << endl
    << "labelloc=t;" << endl;
    
    for (size_t cond=0; cond < getTransTables(0).size(); cond ++) {
        saveDOTCondition(out, cond, colors, showAllNodes, mergeNodes, used, isParent);
    }

    saveDOTNodes(out, used, isParent);
    
	out << endl << "}";
	out.close();
}

void EnvGLNCmp::saveDotFilek2(const string & filename)
{ //support visualization with different parents under two conditions

	ofstream out(filename.c_str());
	out.precision(2);

	//green,blue,chocolate represents common,absolute and relative differential respectively
	// string colors[3] = {"green","blue","chocolate"};

	string colors[] = {"steelblue1","chocolate", "palegreen", "gray"};

	out << "digraph GLNComparativeModeling {" << endl
		<< "ratio=auto;" << endl
		<< "margin=0;" <<endl
        << "fontname=\"Helvetica\";" << endl
        << "edge[fontname=\"Helvetica\"];" << endl
		<< "node[fontname=\"Helvetica\"];" << endl
        << "node[fontcolor=\"black\"];" << endl
		<< "rankdir=TB;" << endl
		<< "labelloc=t;" << endl;
	
	//condition 1
	out << "subgraph \"cluster_condition1\" {" << endl;
	out << "label=\"" << "first condition\";" << endl;
    
    bool showAllNodes = false; // Display all singleton nodes (true) or not (false)
    
    vector<bool> used(getGLNSize(), showAllNodes);
    
	for(size_t i=0; i < getGLNSize(); i++) {
        string nodeNameDest = quote_as_needed(getDiffGLN().getNodeName(i)+"C1");
        if((getType(i)==REL_DIFF||getType(i)==ABS_DIFF)) {
            //first condition
            for(size_t j=0; j<getTransTables(i)[0].getParents().size(); j++)
            {
                int src = getTransTables(i)[0].getParents()[j];
                string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1)+"C1");
                out	<< nodeNameSrc << " -> " << nodeNameDest;
                
                int penwidth = choose_penwidth( getDiffTransTable(i).getpValue() );
                
                out << " [arrowhead=normal,penwidth=" << penwidth
                << ",color="<<colors[0]<<",style=solid, ];" << endl;
                
                used[src-1] = true;  // the node is used for an edge
            }
            
            used[i] = true;
            
        } else if(getType(i)==CONSERVED) {
            
            vector<int> parents1 = getTransTables(i)[0].getParents();
            
            for(size_t j=0; j<parents1.size(); j++)
            {
                int src = parents1[j];
                string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1)+"C1");
                out	<< nodeNameSrc << " -> " << nodeNameDest;
                
                int penwidth = choose_penwidth( getPooledTransTable(i).getpValue() );
                
                out << " [arrowhead=normal,penwidth=" << penwidth
                << ",color="<<colors[2]<<",style=solid, ];" << endl;
                
                used[src-1] = true;  // the node is used for an edge
            }
            
            used[i] = true;
            
        } /*else {  //as discussed with Joe, do not need to display null edges
           //first condition
           for(size_t j=0; j<getTransTables(i)[0].getParents().size(); j++)
           {
           int src = getTransTables(i)[0].getParents()[j];
           string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1)+"C1");
           out	<< nodeNameSrc << " -> " << nodeNameDest;
           
           int penwidth = 1;
           
           out << " [arrowhead=normal,penwidth=" << penwidth
           << ",color="<<colors[3]<<",style=solid, ];" << endl;
           }
           
           }*/
	}
    
	for(size_t i=0; i < getGLNSize(); i++) {
        
        if(used[i]) {
            string nodeNameDest = quote_as_needed(getDiffGLN().getNodeName(i)+"C1");
            out << nodeNameDest;
            out<<  "[label="<< quote_as_needed(getDiffGLN().getNodeName(i)) <<",fontcolor=white,color=transparent,shape=ellipse,style=filled,fillcolor="<<colors[0]<<"];";
            out << endl;
        }
	}
    
	out << "}" << endl << endl;
    
	//condition 2
	out << "subgraph \"cluster_condition2\" {" << endl;
	out << "label=\"" << "second condition\";" << endl;
	for(int i=0; i < (int) getGLNSize(); i++) {
        string nodeNameDest = quote_as_needed(getDiffGLN().getNodeName(i)+"C2");
        if((getType(i)==REL_DIFF||getType(i)==ABS_DIFF)) {
            //first condition
            for(size_t j=0; j<getTransTables(i)[1].getParents().size(); j++)
            {
                int src = getTransTables(i)[1].getParents()[j];
                string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1)+"C2");
                out	<< nodeNameSrc << " -> " << nodeNameDest;
                
                int penwidth = choose_penwidth( getDiffTransTable(i).getpValue() );
                
                out << " [arrowhead=normal,penwidth=" << penwidth
                << ",color="<<colors[1]<<",style=solid, ];" << endl;
            }
            
        } else if(getType(i)==CONSERVED) {
            
            vector<int> parents1 = getTransTables(i)[1].getParents();
            
            for(size_t j=0; j<parents1.size(); j++)
            {
                int src = parents1[j];
                string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1)+"C2");
                out	<< nodeNameSrc << " -> " << nodeNameDest;
                
                int penwidth = choose_penwidth( getPooledTransTable(i).getpValue() );
                
                out << " [arrowhead=normal,penwidth=" << penwidth
                << ",color="<<colors[2]<<",style=solid, ];" << endl;
            }
            
        } /*else {
           //first condition
           for(size_t j=0; j<getTransTables(i)[1].getParents().size(); j++)
           {
           int src = getTransTables(i)[1].getParents()[j];
           string nodeNameSrc = quote_as_needed(getDiffGLN().getNodeName(src-1)+"C2");
           out	<< nodeNameSrc << " -> " << nodeNameDest;
           
           int penwidth = 1;
           
           out << " [arrowhead=normal,penwidth=" << penwidth
           << ",color="<<colors[3]<<",style=solid, ];" << endl;
           }
           
           }*/
	}
	for(int i=0; i < (int) getGLNSize(); i++)
	{	
		string nodeNameDest = quote_as_needed(getDiffGLN().getNodeName(i)+"C2");
		out << nodeNameDest;
		out<<  "[label="<< quote_as_needed(getDiffGLN().getNodeName(i)) <<",fontcolor=white,color=transparent,shape=ellipse,style=filled,fillcolor="<<colors[1]<<"];";
		out << endl;
	}			
	out << "}" << endl << endl;
    
	out << endl << "}";
	out.close();
}

void EnvGLNCmp::saveTopology(const string & fileTopology, int typeIntx, int typeParents, int typeChild) const
{
	if(fileTopology.empty()) return;
    
	// Pathway topology( m_trajCols[0].nNodes(), m_trajCols.size() );
    
    Topology topology( m_trajCols[0].getIntNodeNames(), m_trajCols.size() );
    
	for(int i=0; i < (int) getGLNSize(); i++) {
		
		if( getType(i) == typeIntx 
           && m_childWorkingZones[i] == typeChild 
           && m_workingZones[i] == typeParents ) {
            
            vector<int> parents;
            
            switch(typeIntx) {
				case NULL_INTX:
					continue;
                    
				case CONSERVED:
					parents = m_pooledTransTables[i].getParents();
					break;
                    
				case REL_DIFF:
				case ABS_DIFF:
					parents = m_diffTransTables[i].getParents();
					break;
            }
            topology.set(i, parents);
		}
        
	}
    
	// topology.save( fileTopology, m_trajCols[0].getIntNodeNames() );
	topology.save( fileTopology );
    
}

void EnvGLNCmp::saveTopology(const string & fileTopology) const
{ // Save all non-null interactions where both child and parents have changed
  // working zones
    
	if(fileTopology.empty()) return;
    
    Topology topology( m_trajCols[0].getIntNodeNames(), m_trajCols.size() );
    
	for(int i=0; i < (int) getGLNSize(); i++) {
		
		if( getType(i) != NULL_INTX
           && m_childWorkingZones[i] == CHANGED
           && m_workingZones[i] == CHANGED ) {
            
            vector<int> parents;
            
            switch( getType(i) ) {
				case NULL_INTX:
					continue;
                    
				case CONSERVED:
					parents = m_pooledTransTables[i].getParents();
					break;
                    
				case REL_DIFF:
				case ABS_DIFF:
					// parents = m_diffTransTables[i].getParents();
                    for (size_t k=0; k<m_TransTables[i].size(); k++) {
                        parents.insert(parents.end(),
                                       m_TransTables[i][k].getParents().cbegin(),
                                       m_TransTables[i][k].getParents().cend());
                    }
                    // Remove duplicate parents:
                    sort(parents.begin(), parents.end());
                    vector<int>::iterator it = unique(parents.begin(), parents.end());
                    parents.resize( std::distance(parents.begin(), it) );
					break;
            }
            topology.set(i, parents);
		}
	}
	topology.save( fileTopology );
}

void EnvGLNCmp::saveAllInteractionToDot()
{
    /*
	string pathwayName = m_gc.m_pathway.getFile().size()==0?
    m_gc.m_candidateTopology.getFile():
    m_gc.m_pathway.getFile();
    
	pathwayName = pathwayName.substr(0,pathwayName.find(".txt"));
    
    */
    
	FILE * fp;
    
	// fp = fopen((pathwayName+"DiffPathway.dot").c_str(),"w");
    
	string base;
    if(m_gc.m_pathways.size()==1)
	{
		base = m_gc.m_pathways[0].getFile().substr(0, m_gc.m_pathways[0].getFile().find(".txt"));
	}else if(m_gc.m_candidateTopologys.size()==1){
		base = m_gc.m_candidateTopologys[0].getFile().substr(0, m_gc.m_candidateTopologys[0].getFile().find(".txt"));
	}

    fp = fopen((base + "-Diff.dot").c_str(),"w");
    
	fprintf(fp,"%s\n", "digraph GRN {");
	fprintf(fp,"%s\n", "ratio=auto;");
	fprintf(fp,"%s\n", "margin=0;");
	fprintf(fp,"%s\n", "node[fontname=\"Helvetica\"];");
    fprintf(fp,"%s\n", "fontname=\"Helvetica\";");
	fprintf(fp,"%s\n", "node[fontcolor=\"black\"];");
    
	vector<vector<vector<int> > > pathways;
	if(m_gc.m_pathways.size()==1)
	{
		pathways = m_gc.m_pathways[0].get();
	}else if(m_gc.m_candidateTopologys.size()==1){
		pathways = m_gc.m_candidateTopologys[0].get();
	}
	vector<string> intNodeNames = m_trajCols[0].getIntNodeNames();
	for(size_t i=0; i<pathways.size(); i++)
	{
		for(size_t j=0; j<pathways[i][0].size(); j++)
		{
			//currently only supports 2 experimental conditions
			vector<int> parents1 = getTransTables(i)[0].getParents();
			vector<int> parents2 = getTransTables(i)[1].getParents();
			vector<int> parentIntersection = intersectionSet(parents1,parents2);
			bool flag1 = false; //denote whether node are in parents1
			bool flag2 = false; //denote whether node are in parents2
			bool flag3 = false; //denote whether node are in parents3
			for(size_t m=0; m<parents1.size();m++)
			{
				if(parents1[m]==pathways[i][0][j])
				{
					flag1 = true;
					break;
				}
			}
			for(size_t m=0; m<parents2.size();m++)
			{
				if(parents2[m]==pathways[i][0][j])
				{
					flag2 = true;
					break;
				}
			}
			for(size_t m=0; m<parentIntersection.size();m++)
			{
				if(parentIntersection[m]==pathways[i][0][j])
				{
					flag3 = true;
					break;
				}
			}
            
            string color, style;
            
			if(flag3) {
                color = "brown";
                style = "solid";
                
			} else if(flag1) {        
                color="orange";
                style = "solid";

			} else if(flag2) {
                color="green1";
                style = "solid";
                
			} else {
				color="grey";
                style = "dashed";
			}

            fprintf(fp, "%s\t%s\t%s\t%s%s%s%s%s\n", 
                    quote_as_needed(intNodeNames[pathways[i][0][j]-1]).c_str(),
                    " -> ", 
                    quote_as_needed(intNodeNames[i]).c_str(), 
                    //"[label=-1,
                    "[arrowhead=normal,color=", color.c_str(), ",style=", style.c_str(), " ];");

            /*
			if(flag3)
			{
				fprintf(fp,"%s%s%s\t%s\t%s%s%s\t%s\n", "\"",intNodeNames[pathways[i][0][j]-1].c_str()),"\"",
                        "->","\"",intNodeNames[i].c_str(),"\"",
                        "[label=-1,arrowhead=normal,color=brown,style=solid, penwidth=2];");
			}else if(flag1)
			{
				fprintf(fp,"%s%s%s\t%s\t%s%s%s\t%s\n", "\"",intNodeNames[pathways[i][0][j]-1].c_str(),"\"",
                        "->","\"",intNodeNames[i].c_str(),"\"",
                        "[label=-1,arrowhead=normal,color=orange,style=solid, penwidth=2];");
			}else if(flag2)
			{
				fprintf(fp,"%s%s%s\t%s\t%s%s%s\t%s\n", "\"",intNodeNames[pathways[i][0][j]-1].c_str(),"\"",
                        "->","\"",intNodeNames[i].c_str(),"\"",
                        "[label=-1,arrowhead=normal,color=green1,style=solid, penwidth=2];");
			}else
			{
				fprintf(fp,"%s%s%s\t%s\t%s%s%s\t%s\n", "\"",intNodeNames[pathways[i][0][j]-1].c_str(),"\"",
                        "->","\"",intNodeNames[i].c_str(),"\"",
                        "[label=-1,arrowhead=normal,color=grey,style=dashed ];");
			}
            */
            
		}
	}
    
	//type 0 denotes null interation, 1 denotes relative differential, 2 denotes absolute differential, 3 denotes homogenous interaction
	vector<int> nodesOnPathway(pathways.size(), 0);
    
	for(size_t i=0; i<pathways.size(); i++)
	{
		if(pathways[i][0].size()>0)
		{	
			nodesOnPathway[i] = 1;
			for(size_t j=0;j<pathways[i][0].size();j++)
			{
				nodesOnPathway[pathways[i][0][j]-1] = 1;
			}
		}
        
	}
	vector<double> pvalues;
	for(size_t i=0; i<nodesOnPathway.size(); i++)
	{			
		if(nodesOnPathway[i] == 1)
		{
			pvalues.push_back(getDiffTransTable(i).getpValue());
		}
	}
	sort(pvalues.begin(),pvalues.end());
	for(size_t i=0; i<nodesOnPathway.size(); i++)
	{			
		if(nodesOnPathway[i] == 1)
		{
			double pvalue = getDiffTransTable(i).getpValue();
            string name = quote_as_needed(intNodeNames[i]);
			if(pvalue<pvalues[pvalues.size()/4])
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=indianred1];");
			}else if(pvalue<pvalues[pvalues.size()/2])
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=indianred2];");					
			}else if(pvalue<pvalues[pvalues.size()*3/4])
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=indianred3];");					
			}else
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=indianred4];");					
			}
		}
	}
	fprintf(fp,"%s\n","}");
	fclose(fp);
    
	//the following code is to visualized pathway by conserved interactions
	// fp = fopen((pathwayName+"CommPathway.dot").c_str(),"w");
    
    fp = fopen(( base + "-Conserved.dot" ).c_str(),"w");
    
	fprintf(fp,"%s\n", "digraph GRN {");
	fprintf(fp,"%s\n", "ratio=auto;");
	fprintf(fp,"%s\n", "margin=0;");
	fprintf(fp,"%s\n", "node[fontname=\"Helvetica\"];");
    fprintf(fp,"%s\n", "fontname=\"Helvetica\";");

	fprintf(fp,"%s\n", "node[fontcolor=\"black\"];");
    
	for(size_t i=0; i<pathways.size(); i++)
	{
		for(size_t j=0; j<pathways[i][0].size(); j++)
		{
			vector<int> parents1 = getTransTables(i)[0].getParents();
			vector<int> parents2 = getTransTables(i)[1].getParents();
			vector<int> parentIntersection = intersectionSet(parents1,parents2);
			// bool flag1 = false; //denote whether node are in parents1
			// bool flag2 = false; //denote whether node are in parents2
			bool flag3 = false; //denote whether node are in parents3
			for(size_t m=0; m<parents1.size();m++)
			{
				if(parents1[m]==pathways[i][0][j])
				{
					// flag1 = true;
					break;
				}
			}
			for(size_t m=0; m<parents2.size();m++)
			{
				if(parents2[m]==pathways[i][0][j])
				{
					// flag2 = true;
					break;
				}
			}
			for(size_t m=0; m<parentIntersection.size();m++)
			{
				if(parentIntersection[m]==pathways[i][0][j])
				{
					flag3 = true;
					break;
				}
			}
            
			if(flag3)
			{
				fprintf(fp,"%s\t%s\t%s\t%s\n", quote_as_needed(intNodeNames[pathways[i][0][j]-1]).c_str(), 
                        "->", quote_as_needed(intNodeNames[i]).c_str(),
                        //"[label=-1,
                        "[arrowhead=normal,color=blue,style=solid, penwidth=2];");
			}else
			{
				fprintf(fp,"%s\t%s\t%s\t%s\n", quote_as_needed(intNodeNames[pathways[i][0][j]-1]).c_str(),
                        "->", quote_as_needed(intNodeNames[i]).c_str(),
                        //"[label=-1,
                        "[arrowhead=normal,color=grey,style=dashed];");
			}
		}
	}
    
	pvalues.clear();
	for(size_t i=0; i<nodesOnPathway.size(); i++)
	{			
		if(nodesOnPathway[i] == 1)
		{
			pvalues.push_back(getPooledTransTable(i).getpValue());
		}
	}
	sort(pvalues.begin(),pvalues.end());
	for(size_t i=0; i<nodesOnPathway.size(); i++)
	{			
		if(nodesOnPathway[i] == 1)
		{
			double pvalue = getPooledTransTable(i).getpValue();
            
            string name = quote_as_needed(intNodeNames[i]);
            
			if(pvalue<pvalues[pvalues.size()/4])
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=blue1];");
			}else if(pvalue<pvalues[pvalues.size()/2])
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=blue2];");
			}else if(pvalue<pvalues[pvalues.size()*3/4])
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=blue3];");
			}else
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=blue4];");
			}
		}
	}
	fprintf(fp,"%s\n","}");
	fclose(fp);
    
	//the following code is to visualized pathway by working zone change
	// fp = fopen((pathwayName+"WorkingZonePathway.dot").c_str(),"w");
    
    fp = fopen((base+"-WkZn.dot").c_str(), "w");
    
	fprintf(fp,"%s\n", "digraph GRN {");
	fprintf(fp,"%s\n", "ratio=auto;");
	fprintf(fp,"%s\n", "margin=0;");
	fprintf(fp,"%s\n", "node[fontname=\"Helvetica\"];");
    fprintf(fp,"%s\n", "fontname=\"Helvetica\";");

	fprintf(fp,"%s\n", "node[fontcolor=\"black\"];");
    
	for(size_t i=0; i<pathways.size(); i++)
	{
		for(size_t j=0; j<pathways[i][0].size(); j++)
		{
            
            fprintf(fp,"%s\t%s\t%s\t%s\n", 
                    quote_as_needed(intNodeNames[pathways[i][0][j]-1]).c_str(), 
					"->", 
                    quote_as_needed(intNodeNames[i]).c_str(), 
                    //"[label=-1,
                    "[arrowhead=normal,color=grey,style=dashed];");
		}
	}
    
	pvalues.clear();
	for(size_t i=0; i<nodesOnPathway.size(); i++)
	{			
		if(nodesOnPathway[i] == 1)
		{
			pvalues.push_back(getChildpChisq_z(i));
		}
	}
	sort(pvalues.begin(),pvalues.end());
	for(size_t i=0; i<nodesOnPathway.size(); i++)
	{			
		if(nodesOnPathway[i] == 1)
		{
			double pvalue = getChildpChisq_z(i);
            
            string name = quote_as_needed(intNodeNames[i]);
            
			if(pvalue<pvalues[pvalues.size()/4])
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=wheat1];");
			}else if(pvalue<pvalues[pvalues.size()/2])
			{
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=wheat2];");					
			}else if(pvalue<pvalues[pvalues.size()*3/4]) {
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=wheat3];");					
			} else {
				fprintf(fp,"%s\t%s%3.2e%s\n", name.c_str(),
                        "[shape=ellipse,style=filled,label=\"", pvalue, "\",fillcolor=wheat4];");					
			}
		}
	}
	fprintf(fp,"%s\n","}");
	fclose(fp);
    
}

void parseResultFile(const string & resultFileLocation, vector<double> & Chisq, vector<int> & Df, vector<double> & Pvalue)
//added by Yang Zhang 7.30.2010 used in Main.cpp for permutation study
// Relocated to EnvGLNCmpIO.cpp by Joe Song
{
    
    //open file for load data check if file is open
    ifstream ifs(resultFileLocation.c_str(), ios::in);
    if(!ifs.is_open()) 
    {
        cerr << "ERROR: openning file \"" << resultFileLocation << "\"" << endl;
        GLNExit(EXIT_FAILURE);
    }
    
    string line; 
    
	size_t index = 0;

	while(getline(ifs,line))
	{
		line = line.substr(line.find_first_of('\t')+1);
		Pvalue[index] = atof(line.substr(0,line.find_first_of('\t')).c_str());
		
		line = line.substr(line.find_first_of('\t')+1);
		Df[index] = atoi(line.substr(0,line.find_first_of('\t')).c_str());

		line = line.substr(line.find_first_of('\t')+1);
		Chisq[index] = atof(line.substr(0,line.find_first_of('\t')).c_str());
		
		index++;
	}
/*
    getline(ifs,line); //Chi-d header
    getline(ifs,line); //p_d
    Pvalue[0] = atof(line.substr(line.find('')+1).c_str());
    getline(ifs,line); //df_d
    Df[0] = atoi(line.substr(line.find(':')+1).c_str());
    getline(ifs,line); //chisq_d
    Chisq[0] = atof(line.substr(line.find(':')+1).c_str());
    
    getline(ifs,line); //Chi-c header
    getline(ifs,line); //p_c
    Pvalue[1] = atof(line.substr(line.find(':')+1).c_str());
    getline(ifs,line); //df_c
    Df[1] = atoi(line.substr(line.find(':')+1).c_str());
    getline(ifs,line); //chisq_c
    Chisq[1] = atof(line.substr(line.find(':')+1).c_str());
    
    getline(ifs,line); //Chi_t header
    getline(ifs,line); //p_t
    Pvalue[2] = atof(line.substr(line.find(':')+1).c_str());
    getline(ifs,line); //df_t
    Df[2] = atoi(line.substr(line.find(':')+1).c_str());
    getline(ifs,line); //chisq_t
    Chisq[2] = atof(line.substr(line.find(':')+1).c_str());
    
    getline(ifs,line); //Chi_parent header
    getline(ifs,line); //p_p_z
    Pvalue[3] = atof(line.substr(line.find(':')+1).c_str());
    getline(ifs,line); //df_p_z
    Df[3] = atoi(line.substr(line.find(':')+1).c_str());
    getline(ifs,line); //chisq_p_z
    Chisq[3] = atof(line.substr(line.find(':')+1).c_str());
    
    getline(ifs,line); //Chi_child header
    getline(ifs,line); //p_c_z
    Pvalue[4] = atof(line.substr(line.find(':')+1).c_str());
    getline(ifs,line); //df_c_z
    Df[4] = atoi(line.substr(line.find(':')+1).c_str());
    getline(ifs,line); //chisq_c_z
    Chisq[4] = atof(line.substr(line.find(':')+1).c_str());
  */  
}

void EnvGLNCmp::recordResults(int child, const TransitionTable &tt, string description)
// MS. Dec 4, 2011 const vector<Node> &nodes, 
{ 
	//m_resultrecord.push_back(_result);
	stringstream out1;  //child node id
	//convert int to string
	out1 << child + 1;
	string _result = description + " " + out1.str();
    //string _result = out1.str();
    
	_result = _result + " " + getNodeName(child);
	_result = _result + " " + getNodeType(child);
    
	stringstream out2; //parent size
	out2 << tt.getParents().size();
	_result = _result + " " + out2.str() + " ";
    
	for(size_t j=0; j < tt.getParents().size(); j++) 
	{
		//convert int to string
		stringstream out3; //parent node id
		out3 << tt.getParents()[j];
		_result = _result + out3.str() + ":" + getNodeName(tt.getParents()[j]-1) + ",";
		//_result = _result + out3.str() +  ",";
	}
    
	stringstream out4;
	out4<< tt.getpValue();
	_result = _result + " " + out4.str();
    
	stringstream out5;
	out5<< tt.getChisq();
	_result = _result + " " + out5.str();
    
	stringstream out6;
	out6<<tt.getDf();
	_result = _result + " " + out6.str();
    
	m_recordresult[child].push_back(_result.c_str());
	m_recordcounts[child].push_back(tt.countsToString());
}


void EnvGLNCmp::recordResults(int child, const TransitionTable &ttHet, const TransitionTable &ttHom, const TransitionTable &ttTot, string description)
// MS. Dec 4, 2011 const vector<Node> &nodes, 
{
	//m_resultrecord.push_back(_result);
	stringstream out1;  //child node id
	//convert int to string
	out1 << child + 1;
	string _result = description + " " + out1.str();
	//string _result = out1.str();

	_result = _result + " " + getNodeName(child);
	_result = _result + " " + getNodeType(child);

	stringstream out2; //parent size
	out2 << ttHet.getParents().size();
	_result = _result + " " + out2.str() + " ";

	for (size_t j = 0; j < ttHet.getParents().size(); j++)
	{
		//convert int to string
		stringstream out3; //parent node id
		out3 << ttHet.getParents()[j];
		_result = _result + out3.str() + ":" + getNodeName(ttHet.getParents()[j] - 1) + ",";
		//_result = _result + out3.str() +  ",";
	}

	stringstream out4;
	out4 << ttHet.getpValue();
	_result = _result + " " + out4.str();

	stringstream out5;
	out5 << ttHet.getChisq();
	_result = _result + " " + out5.str();

	stringstream out6;
	out6 << ttHet.getDf();
	_result = _result + " " + out6.str();

	stringstream out7;
	out7 << ttHom.getpValue();
	_result = _result + " " + out7.str();

	stringstream out8;
	out8 << ttHom.getChisq();
	_result = _result + " " + out8.str();

	stringstream out9;
	out9 << ttHom.getDf();
	_result = _result + " " + out9.str();

	stringstream out10;
	out10 << ttTot.getpValue();
	_result = _result + " " + out10.str();

	stringstream out11;
	out11 << ttTot.getChisq();
	_result = _result + " " + out11.str();

	stringstream out12;
	out12 << ttTot.getDf();
	_result = _result + " " + out12.str();

	m_recordresult[child].push_back(_result.c_str());
	m_recordcounts[child].push_back(ttHet.countsToString());
}
