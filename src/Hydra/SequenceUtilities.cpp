// ***************************************************************************
// CSequenceUtilities - handles basic sequence manipulation.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "SequenceUtilities.h"

// Performs an in-place reverse complement conversion
void CSequenceUtilities::GetReverseComplement(string &seqBases) {

	// reverse the sequence
	ReverseSequence(seqBases);

	// swap the bases
	for(unsigned int i = 0; i < seqBases.size(); i++) {
		switch(seqBases[i]) {
			case 'A':
				seqBases[i] = 'T';
				break;
			case 'C':
				seqBases[i] = 'G';
				break;
			case 'G':
				seqBases[i] = 'C';
				break;
			case 'T':
				seqBases[i] = 'A';
				break;
			default:
				break;
		}
	}
}
