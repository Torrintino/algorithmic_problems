#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#define RESET "\033[0m"
#define RED "\033[31m"
#define GREEN "\033[32m"

using Sol = std::vector<std::vector<double>>;

struct Bin {
    Bin() = default;

    void addItem(double item) {
        
	    items.push_back(item);
	    weight += item;
	    
	    if(weight > 1) throw std::runtime_error("error - Item can not be added, maximum weight would be exceeded otherwise.");
    }

    double getWeight() const noexcept { return weight; }

    const std::vector<double> &getItems() const noexcept { return items; }

    void clear() {
        items.clear();
        weight = 0;
    }

private:
    std::vector<double> items;
    double weight = 0;
};

Bin &firstAvailableBin(std::vector<Bin> &bins, double item) {
	
	for(int i=0; i<bins.size(); i++)
	    if((bins[i].getWeight() + item) <= 1)
		    return bins[i];

	bins.push_back(Bin());
	return bins.back();
}

std::vector<Bin> firstFit(const std::vector<double> &items) {
	
	std::vector<Bin> res (1);

	for(double item : items)
		firstAvailableBin(res, item).addItem(item);

	return res;
}

std::vector<Bin> firstFitDecreasing(const std::vector<double> &items) {
	
	std::vector<double> toSort = items;

	std::sort(toSort.begin(), toSort.end(), std::greater_equal<double>());

	return firstFit(toSort);
}

// No more tasks below this line

void printBins(const std::vector<Bin> &bins) {
    if (bins.empty())
        std::cout << "No bins\n";

    std::cout << '\n';
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    int rows = std::max_element(
                   bins.begin(), bins.end(),
                   [](auto &b1, auto &b2) { return b1.getItems().size() < b2.getItems().size(); })
                   ->getItems()
                   .size();

    for (int row = rows - 1; row >= 0; --row) {
        for (const Bin &bin : bins) {
            const auto &items = bin.getItems();
            if (items.size() > row)
                std::cout << "| " << items[row] << " |  ";
            else
                std::cout << "          ";
        }
        std::cout << '\n';
    }

    for (const Bin &bin : bins)
        std::cout << "--------  ";
    std::cout << '\n';

    for (const Bin &bin : bins)
        std::cout << "  " << bin.getWeight() << "    ";
    std::cout << "\n\n";
}

class Test {
    static void printFailed() { std::cout << RED << "Failed\n" << RESET; }

    static void printPassed() { std::cout << GREEN << "Passed\n" << RESET; }

public:
    static void doTest(std::vector<double> items, Sol solFF, Sol solFFD) {
        std::cout << "\nTesting [";
        for (size_t i = 0; i < items.size(); ++i) {
            std::cout << items[i];
            if (i != items.size() - 1)
                std::cout << ", ";
        }
        std::cout << "]\n";

        std::cout << "Test add item: ";
        if (!testAddItem(items)) {
            printFailed();
            std::cout << "Stopping current test suite here\n";
            return;
        }
        printPassed();

        std::cout << "Test FF: ";
        auto output = firstFit(items);
        if (testAlgo(output, solFF)) {
            printPassed();
            printBins(output);
        } else {
            printFailed();
        }

        std::cout << "Test FFD: ";
        output = firstFitDecreasing(items);
        if (testAlgo(output, solFFD)) {
            printPassed();
            printBins(output);
        } else {
            printFailed();
        }
    }

private:
    static bool testAlgo(const std::vector<Bin> &output, const Sol &sol) {
        if (output.size() != sol.size())
            return false;

        for (size_t i = 0; i < sol.size(); ++i) {
            auto items = output[i].getItems();
            if (items.size() != sol[i].size())
                return false;

            for (size_t j = 0; j < sol[i].size(); ++j) {
                if (items[j] != sol[i][j])
                    return false;
            }
        }
        return true;
    }

    static bool testAddItem(const std::vector<double> &items) {
        Bin bin{};
        double weight = 0;
        int i = 0, base = 0;

        while (i < items.size()) {
            try {
                bin.addItem(items[i]);
                if (weight + items[i] > 1)
                    return false;
                weight += items[i];
            } catch (const std::runtime_error &err) {
                if (weight + items[i] <= 1)
                    return false;
                bin.clear();
                weight = 0;
                base = i;
                continue;
            }

            for (int j = base; j <= i; ++j) {
                if (items[j] != bin.getItems()[j - base])
                    return false;
            }

            if (weight != bin.getWeight())
                return false;
            ++i;
        }

        return true;
    }
};

int main(int arcg, char **argv) {
    Test::doTest({0.1, 0.2, 0.3, 0.5, 0.1, 0.6, 0.4, 0.4},
                 {{0.1, 0.2, 0.3, 0.1}, {0.5, 0.4}, {0.6, 0.4}},
                 {{0.6, 0.4}, {0.5, 0.4, 0.1}, {0.3, 0.2, 0.1}});
    Test::doTest({0.5, 0.6, 0.7, 0.4, 0.05, 0.06, 0.17, 0.55, 0.23, 0.44, 0.02},
                 {{0.5, 0.4, 0.05, 0.02}, {0.6, 0.06, 0.17}, {0.7, 0.23}, {0.55, 0.44}},
                 {{0.7, 0.23, 0.06}, {0.6, 0.4}, {0.55, 0.44}, {0.5, 0.17, 0.05, 0.02}});
    Test::doTest(
        {0.01, 0.84, 0.16, 0.73, 0.37, 0.08, 0.09, 0.21, 0.95, 0.04, 0.11, 0.63},
        {{0.01, 0.84, 0.08, 0.04}, {0.16, 0.73, 0.09}, {0.37, 0.21, 0.11}, {0.95}, {0.63}},
        {{0.95, 0.04, 0.01}, {0.84, 0.16}, {0.73, 0.21}, {0.63, 0.37}, {0.11, 0.09, 0.08}});
    return 0;
}
