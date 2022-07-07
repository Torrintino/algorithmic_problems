#include <cassert>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "bin-packing.cpp"

class HillClimbing {
    // Current solution
    std::vector<Bin> bins;

    // Vectors used to partition the current solution
    std::vector<Bin> pi, rho;

    // Distribution and random number generator used to generate random numbers between 0 and 1
    std::uniform_real_distribution<double> dist;
    std::mt19937_64 gen;

    // Number of bins used by FF on each step of the hill climbing
    std::vector<int> binSizes;

    // Iterates over all the pairs of items in the given bin. For each pair
    // yields the index and the weight of both items.
    template <class Lambda>
    void forEachItemsPair(const Bin &bin, Lambda lambda) {
        auto &items = bin.getItems();
        for (int i = 0; i < items.size(); ++i)
            for (int j = i + 1; j < items.size(); ++j)
                lambda(i, items[i], j, items[j]);
    }

    // Iterates over all the items of the given bin and yields their index and
    // weight
    template <class Lambda>
    void forEachItem(const Bin &bin, Lambda lambda) {
        auto &items = bin.getItems();
        for (int i = 0; i < items.size(); ++i)
            lambda(i, items[i]);
    }

    // Copies the bins in the 'bins' vector into 'pi' and 'rho'. A bin is
    // copied to 'rho' with probability 1/bins.size(), or to 'pi' otherwise.
    //
    // Note: before starting to add new items to pi and rho, remember to empty
    // them.
    void divide() {

	    double size = static_cast<double>(bins.size());
	    double prob = 1/size;
      	   
	    do
	    {
		    pi.clear();
		    rho.clear();

		    for(Bin currentBin : bins)
		    {
			    double rnd = dist(gen); 
			    

			    if(rnd <= prob) rho.push_back(currentBin);
			    else pi.push_back(currentBin);
		    }

	    }while( rho.empty() || pi.empty() );

    }

    // Runs the improvement phase. Returns true if any bin has been modified,
    // false otherwise.
    bool improve() {

	bool notUpdated = true;	

	for(Bin &binPi: pi)
	{

	//PAIR VS PAIR
	//------------	
		
		forEachItemsPair(binPi, [&](int i, double item_i, int j, double item_j) {
		
				for(Bin &binRho : rho)
				{		
					//FOR EACH ITEM PAIR IN RHO
					forEachItemsPair(binRho, [&](int k, double item_k, int l, double item_l) {

							double sum = item_k + item_l - (item_i + item_j);
							
							if( sum > 0 && (binPi.getWeight() + sum) < 1.0 ) 
							{
								binRho.update(k, item_i);
								binRho.update(l, item_j);

								binPi.update(i, item_k);
								binPi.update(j, item_l);
							
								notUpdated = false;
							}


					});

					
				}	
				
				
				
		});

	//------------
	//PAIR VS PAIR
	



	//PAIR VS SINGLE
	//--------------
	
		forEachItemsPair(binPi, [&](int i, double item_i, int j, double item_j) {
				
				for(Bin &binRho : rho)
				{

					if( j >= binPi.getItems().size() ) break;

					forEachItem(binRho, [&](int k, double item_k) {

							double sum = item_k - (item_i + item_j);
							
							if( sum > 0 && (binPi.getWeight() + sum) < 1.0 )
							{
								binRho.update(k, item_i);
								binRho.addItem(item_j);
								
								binPi.deleteAt(j);
								binPi.update(i, item_k);
									
								notUpdated = false;
							}
							
					});


				}

		});		
		                                                                               
	//--------------
	//PAIR VS SINGLE
	
	//SINGLE vs SINGLE
	//----------------
	
		forEachItem(binPi, [&](int i, double item_i) {
				
				for(Bin &binRho : rho)
				{
					forEachItem(binRho, [&](int k, double item_k) {
							
							double sum = item_k - item_i;
							
							if( sum > 0 && (binPi.getWeight() + sum) < 1.0 )
							{
								binRho.update(k, item_i);

								binPi.update(i, item_k);

								notUpdated = false;
							}			

					});

				}

		});

	//----------------
	//SINGLE VS SINGLE
	
	}


	    return notUpdated;
    }

    // No further tasks below this line.

    void checkDivide() const {
        assert(!pi.empty() && "After divide() pi should not be empty.");
        assert(!rho.empty() && "After divide() rho should not be empty.");
        assert(rho.size() + pi.size() == bins.size());
    }

    void checkWeights() const {
        for (auto &bin : pi)
            assert(bin.getWeight() <= 1 && "Bin weight must not exceed 1");
        for (auto &bin : rho)
            assert(bin.getWeight() <= 1 && "Bin weight must not exceed 1");
    }

public:
    // Constructs an instance of HillClimbing and computes the initial solution using FF.
    HillClimbing(const std::vector<double> &items, int randomSeed = 42) : dist(0, 1) {
        bins = BinPacking::firstFit(items);
        gen.seed(randomSeed);
    }

    // Runs the hill-climbing algorithm, keeps iterating until no bins are changed.
    void climb() {
        bool stop = false;
        do {
            // Divide step
            divide();
            checkDivide();

            // Improve step
            stop = improve();
            checkWeights();

            // Concatenate pi and rho
            std::vector<double> newItems;
            for (Bin &bin : pi)
                newItems.insert(newItems.end(), bin.getItems().begin(), bin.getItems().end());
            for (Bin &bin : rho)
                newItems.insert(newItems.end(), bin.getItems().begin(), bin.getItems().end());

            // Compute new solution
            bins = BinPacking::firstFit(newItems);

            binSizes.push_back(bins.size());
        } while (!stop);
    }

    const std::vector<int> &getBinSizes() const { return binSizes; }
};

// Reads a list of items from a binary file
std::vector<double> readItems(const std::string &path) {
    std::ifstream file(path, std::ios::binary);
    std::vector<double> vec;
    double val;
    while (file.read(reinterpret_cast<char *>(&val), sizeof(double)))
        vec.push_back(static_cast<double>(val));
    return vec;
}

// Stores intermediate solutions (number of bins) in a binary file
template <class T>
void vectorToFile(const std::vector<T> &vec, const std::string &path) {
    std::ofstream outStream(path, std::ios::out | std::ios::binary);
    outStream.write(reinterpret_cast<const char *>(&vec[0]), vec.size() * sizeof(T));
    outStream.close();
}

int main(int arcg, char **argv) {
    for (int i = 1; i <= 5; ++i) {
        auto items = readItems("./data/items_" + std::to_string(i) + ".txt");
        HillClimbing hc(items);
        hc.climb();
        vectorToFile(hc.getBinSizes(), "./data/bin_sizes_" + std::to_string(i) + ".txt");
    }
    return 0;
}
