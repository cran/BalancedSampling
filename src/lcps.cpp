#include <Rcpp.h>
#include <math.h>

//**********************************************
// Author: Wilmer Prentius
// Last edit: 2022-06-28
// Licence: GPL (>=2)
//**********************************************

namespace WP {
typedef unsigned int uint;
}

int intN(WP::uint n) {
  return std::floor(R::runif(0, 1) * n);
}

struct Object {
  WP::uint index;
  WP::uint rindex;
  double probability;
  short int selected;
  double random;
  Object() : index(0), rindex(0), probability(0.0), selected(-1), random(0.0) {}
  // Object(WP::uint idx) : index(idx), rindex(idx), probability(0.0), selected(-1) {}
  
  double selectObjectByRandom() {
    if (random <= probability) {
      selected = 1;
      return 1.0;
    }
    
    selected = 0;
    return 0.0;
  }
  
  void updateSelection(const double prob) {
    probability = prob;
    
    if (prob - 1e-9 <= 0.0) {
      selected = 0;
      return;
    }
    
    if (prob + 1e-9 >= 1.0) {
      selected = 1;
      return;
    }
    
    return;
  }
  
  double updateProbability(Object *o1, const double weight, const double availWeight, const double slag) {
    // What happens to pj?
    // Ii            |  0  |    1
    // --------------|-----|-----
    // pi + pj <= 1  |  w  |    0
    // pi + pj >= 1  |  1  |  1-w
    
    // availWeight is the maximum weight left to allocate,
    // i.e. it is either the remaining weight
    // or the remaining weight divided by the number of equidistant objects
    // weight is the maxweight possible for a unit
    
    // Updating should be according to the smallest of these two weight
    // If the smallest is the maxweight, then the table above is true, and we can
    // simplify some calculations (albeit the if's might not make it simpler in reality :)
    // However, if we're constrained by the amount of weight left to allocate,
    // the regular assignment must be used.
    
    // The rationale behind these if's is to reduce the rounding errors,
    // however it is not tested whether or not this actually has an effect.
    double retWeight;
    if (weight <= availWeight) {
      retWeight = weight;
      
      if (o1->probability + probability <= 1.0) {
        // First row
        if (o1->selected == 1) {
          updateSelection(0.0);
        } else {
          updateSelection(weight);
        }
      } else {
        // Second row
        if (o1->selected == 1) {
          updateSelection(1.0 - weight);
        } else {
          updateSelection(1.0);
        }
      }
      
      return retWeight;
    }
    
    retWeight = availWeight;
    updateSelection(probability - slag * availWeight);
    return retWeight;
  }
};

struct DiPair {
  double distance;
  double weight;
  WP::uint index;
  DiPair() : distance(0.0), weight(0.0), index(0) {}
  void set(double d, double w, WP::uint i) {
    distance = d;
    weight = w;
    index = i;
    return;
  }
};

struct Experiment {
  Rcpp::NumericMatrix &data;
  WP::uint length;
  
  std::vector<Object> &objects;
  std::vector<WP::uint> &indices;
  
  WP::uint unresolvedObjects;
  std::vector<DiPair> &temp_di1;
  std::vector<DiPair> &temp_di2;
  
  Experiment(
    Rcpp::NumericMatrix &dt,
    WP::uint len,
    std::vector<Object> &objs,
    std::vector<WP::uint> &idx,
    std::vector<DiPair> &di1,
    std::vector<DiPair> &di2
  ) :
    data(dt),
    length(len),
    objects(objs),
    indices(idx),
    unresolvedObjects(len),
    temp_di1(di1),
    temp_di2(di2)
  {}
  
  double euclideanDistance(const WP::uint a, const WP::uint b) {
    double d = 0.0;
    double t;
    WP::uint cols = data.ncol();
    
    for (WP::uint j = 0; j < cols; j++) {
      t = data(a, j) - data(b, j);
      d += t * t;
    }
    
    return d;
  }
  
  void setSortTemp_di(const WP::uint current, const WP::uint id) {
    WP::uint tempid;
    WP::uint unresi  = 0;
    
    for (WP::uint i = 0; i < unresolvedObjects; i++) {
      if (current == i) { continue; }
      tempid = indices[i];
      temp_di1[unresi].set(euclideanDistance(id, tempid), 0.0, tempid);
      unresi++;
    }
    
    if (unresi > 1) {
      std::sort(
        temp_di1.begin(),
        temp_di1.begin()+unresi,
        [](DiPair &a, DiPair &b) -> bool {return a.distance < b.distance;}
      );
    }
  }
  
  double maxweight(const WP::uint a, const WP::uint b) {
    Object &o1 = objects[a];
    Object &o2 = objects[b];
    if (o1.probability <= 0.0 || o1.probability >= 1.0) { return 0.0; }
    
    if (o1.probability + o2.probability <= 1.0) {
      return o2.probability / (1.0 - o1.probability);
    }
    
    return (1.0 - o2.probability) / o1.probability;
  }
  
  double distance(const WP::uint current) {
    WP::uint unres = unresolvedObjects - 1;
    WP::uint id = indices[current];
    
    // calcualte distances between objects
    // use temp_di
    setSortTemp_di(current, id);
    
    double oneenough = 1.0 - 1e-9;
    double weight = 0.0;
    
    for (WP::uint i = 0; i < unres; i++) {
      weight += maxweight(id, temp_di1[i].index);
      
      if (weight > oneenough) {
        return temp_di1[i].distance;
      }
    }
    
    return temp_di1[unres-1].distance;
  }
  
  WP::uint chooseObject() {
    if (unresolvedObjects <= 2) {
      return WP::uint(0);
    }
    
    // Get the distance of affection for each object
    for (WP::uint i = 0; i < unresolvedObjects; i++) {
      temp_di2[i].set(distance(i), 0.0, i);
    }
    
    std::sort(
      temp_di2.begin(),
      temp_di2.begin() + unresolvedObjects,
      [](DiPair &a, DiPair &b) -> bool {return a.distance < b.distance;}
    );
    
    WP::uint n;
    double md = temp_di2[0].distance;
    for (n = 1; n < unresolvedObjects; n++) {
      if (temp_di2[n].distance > md) { break; }
    }
    
    if (n == 1) { return temp_di2[0].index; }
    
    return temp_di2[intN(n)].index;
  }
  
  void decideObject(const WP::uint id, const WP::uint rid) {
    unresolvedObjects--;
    WP::uint temp = indices[unresolvedObjects];
    indices[unresolvedObjects] = indices[rid];
    indices[rid] = temp;
    objects[id].rindex = unresolvedObjects;
    objects[temp].rindex = rid;
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector lcps(Rcpp::NumericVector &prob, Rcpp::NumericMatrix &x) {
  WP::uint rows = WP::uint(x.nrow());
  if (prob.length() != rows) {
    Rcpp::stop("length of prob must be same as number of rows in x");
  }
  
  std::vector<Object> objects(rows);
  std::vector<WP::uint> indices(rows);
  std::vector<DiPair> temp_di1(rows);
  std::vector<DiPair> temp_di2(rows);
  Experiment lcps(
      x,
      rows,
      objects,
      indices,
      temp_di1,
      temp_di2
  );
  
  for (WP::uint ii = 0; ii < rows; ii++) {
    if (prob(ii) <= 0.0 || 1.0 <= prob(ii)) {
      Rcpp::stop("all elements of prob must be in (0, 1)");
    }
    
    objects[ii].rindex = ii;
    objects[ii].index = ii;
    objects[ii].random = R::runif(0, 1);
    objects[ii].probability = prob(ii);
    indices[ii] = ii;
  }
  
  WP::uint id, rid, temp_di_len, equals;
  Object *o1, *o2;
  double slag, totalWeight, weight;
  
  // Run loop while there are unresolved objects
  while (lcps.unresolvedObjects > 1) {
    rid = lcps.chooseObject();
    id = lcps.indices[rid];
    o1 = &lcps.objects[id];
    
    // decide if object should be selected
    // slag is set to the selection status of the object
    slag = o1->selectObjectByRandom();
    
    // swaps the places in e.indicies
    lcps.decideObject(id, rid);
    
    // might be unnecessary this check? Not if N=1? But then this should be unreachable?
    if (lcps.unresolvedObjects == 0) {
      break;
    }
    
    // Now slag is set to 1-pi (see eq. 3 in Spatially balanced Poisson sampling, Grafström 2012)
    // i.e. the probability that must be offset from other objects
    slag -= o1->probability;
    // Get a sorted temp_di of objects, i.e. a list of object-indices ordered by distance to o1
    lcps.setSortTemp_di(lcps.unresolvedObjects, id);
    
    // Loop through temp_di and decide objects
    totalWeight = 1.0; // reset totalWeight
    temp_di_len = lcps.unresolvedObjects; // needs to be set as uO changes during loop, but temp_di is "fixed" size
    for (WP::uint i = 0; i < temp_di_len && totalWeight > 1e-9;) {
      // set maxweight of current unit
      lcps.temp_di1[i].weight = lcps.maxweight(id, lcps.temp_di1[i].index);
      
      // loop through possible units w/ equal distance
      for (equals = 1; equals + i < temp_di_len; equals++) {
        // Allow for some floaterror
        if (lcps.temp_di1[i + equals].distance - 1e-9 > lcps.temp_di1[i + equals - 1].distance) { break; }
        
        // set maxweight for any unit w/ equal distance
        lcps.temp_di1[i + equals].weight = lcps.maxweight(id, lcps.temp_di1[i + equals].index);
      }
      
      // if needed, sort the temp_di w/ equal distance so that low weight go before high
      // since [i, i+equals) should only contain equal distance,
      // it should be sufficient to sort according to weight
      if (equals > 1) {
        std::sort(
          lcps.temp_di1.begin() + i,
          lcps.temp_di1.begin() + i + equals,
          [](DiPair &a, DiPair &b) -> bool {return a.weight < b.weight;}
        );
      }
      
      // loop through all elements of equal distance
      for (WP::uint j = i; j < i + equals; j++) {
        o2 = &lcps.objects[lcps.temp_di1[j].index];
        
        weight = o2->updateProbability(
          o1,
          lcps.temp_di1[j].weight,
          totalWeight / double(i + equals - j),
          slag
        );
        
        totalWeight -= weight;
        
        if (o2->selected != -1) {
          lcps.decideObject(lcps.temp_di1[j].index, o2->rindex);
        }
      }
      
      i += equals;
    }
    
    // if we're left with only one object, we should end it here
    if (lcps.unresolvedObjects == 1) {
      id = lcps.indices[0];
      lcps.objects[id].selectObjectByRandom();
    }
  }
  
  int n = int(ceil(Rcpp::sum(prob)));
  std::vector<int> s(n);
  WP::uint count = 0;
  
  for (WP::uint i = 0; i < rows; i++) {
    if (lcps.objects[i].selected != 1) { continue; }
    s[count] = int(i + 1);
    count++;
  }
  
  s.resize(count);
  
  return Rcpp::wrap(s);
}
