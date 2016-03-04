#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hgn.h"
#include <inttypes.h>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <iterator>

struct MinHalo {
  PartId haloid;                /*!< The haloes ID */
  PartId partcount;     /*!< The number of particles in the halo */
  PartDesc *partlist;     /*!< The allocated list of particles  */
  MinHalo(HGNHalo *h) : haloid(h->haloid), partcount(h->partcount),
	partlist(h->partlist) {};
};
typedef	std::vector<MinHalo *> MHL;

/* basic sort - on size of particles*/
static bool sortfnc(MinHalo *h1, MinHalo *h2) 
{
	return h1->partcount < h2->partcount;
}
std::set<PartId> dupids;

std::ostream &operator<<(std::ostream &os, const PartDesc &p) {
    os << p.id << " " << p.mass << " " << p.age << " " << p.z;
    return os;
}

struct PDescLess {
    bool operator()(const PartDesc &p1, const PartDesc &p2) {
	return p1.id < p2.id;
    }
};


/* print a halo discarding duplicates */
void printhalo(MinHalo *mh)
{
	std::vector<PartDesc> ids;

	ids.reserve(mh->partcount); // we'll need this amount of space 
	for(size_t i = 0; i < mh-> partcount; i++) {
		// if insertion .second is true, its a new entry
	    if (dupids.insert(mh->partlist[i].id).second == true)
		ids.push_back(mh->partlist[i]);
	}
	// sort them - might help with locality of reference
	std::sort(ids.begin(), ids.end(), PDescLess());
	// halo header
  //std::cout << "#" << std::endl;
	std::cout << ids.size() << " " << mh->haloid << std::endl;
	// stream out the list
	std::copy(ids.begin(), ids.end(), 
		  std::ostream_iterator<PartDesc>(std::cout, "\n"));
}


/*! remove duplicate ids from haloes */
int main(int argc, char **argv)
{
	HGN *hgn;
	HGNHalo *halo;

	if (argc < 2) {
		fprintf(stderr, "Usage: %s filename|-\n", argv[0]);
		exit(1);
	}
	if (strcmp(argv[1], "-") == 0)
		hgn = openhaloidstream(stdin);
	else
		hgn = openhaloidfile(argv[1]);
	if (hgn == NULL) {
		fprintf(stderr, "Can't open ID file %s\n", argv[1]);
		exit(1);
	}
	/* gather all the haloes */
	MHL halolist;
	while ((halo = getnexthalo(hgn)) != NULL) {
		MinHalo *mh = new MinHalo(halo);
		free(halo);
		std::clog << "Read " << mh->haloid << std::endl;
		halolist.push_back(mh);
	}
	/* sort them */
  std::sort(halolist.begin(), halolist.end(), sortfnc);

	/* output them */
	MHL::iterator iter;
	// header
//	std::cout << "# " << hgn->idline << std::endl;
//	std::cout << "# level: " << hgn->level << std::endl;
//	std::cout << "# author: " << hgn->author << std::endl;
//	std::cout << "# halo finder: "<<hgn->hf << std::endl;
//	if (hgn->date)
//		std::cout << "# date: "<<hgn->date << std::endl;

	// data
	size_t n = 0;
  std::cout << halolist.size() << std::endl;
	for(iter = halolist.begin(); iter != halolist.end(); ++iter) {
		std::clog << n++ << " of " << halolist.size() <<  std::endl;
		printhalo(*iter);
		free((*iter)->partlist);
		(*iter)->partlist = 0;
	}
	exit(0);
}
