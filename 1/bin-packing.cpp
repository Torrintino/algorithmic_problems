#include <stdexcept>
#include <vector>

struct Bin {
    Bin() = default;

    void addItem(double item) {
        if (weight + item > 1)
            throw std::runtime_error("Error, bin weight would exceed 1");
        items.push_back(item);
        weight += item;
    }

    // Updates the item at index `index` to `newItem`.
    // Updates the `weight` variable accordingly.
    void update(int index, double newItem) {

	weight+= (newItem - items[index]);
	items[index] = newItem;	
    }

    // Deletes the item at index `index`.
    // Updates the `weight` variable accordingly.
    void deleteAt(int index) {

        weight-= items[index];
	items.erase(items.begin() + index);
    }

    // No further tasks below this line.

    double getWeight() const noexcept { return weight; }

    const std::vector<double> &getItems() const noexcept { return items; }

    bool empty() const noexcept { return items.empty(); }

    double at(int index) const { return items.at(index); }

    void clear() {
        items.clear();
        weight = 0;
    }

private:
    std::vector<double> items;
    double weight = 0;
};

class BinPacking {

    static Bin &firstAvailableBin(std::vector<Bin> &bins, double item) {
        for (Bin &bin : bins)
            if (bin.getWeight() + item <= 1)
                return bin;

        bins.emplace_back();
        return bins.back();
    }

public:
    static std::vector<Bin> firstFit(const std::vector<double> &items) {
        std::vector<Bin> bins;
        for (double item : items) {
            auto &bin = firstAvailableBin(bins, item);
            bin.addItem(item);
        }

        return bins;
    }
};
