#include <string>
#include <map>
#include <vector>
using namespace std;

struct kmer_models {
    string kmers;
    float mean;
    float stdv;
};

struct kmers_parse {
    string kmers;
    vector<float> samples;
    string contig;
    string name;
};

struct kmers_parse_r {
    vector<pair<int, kmers_parse>> kmers_lines;
    int counter;
    vector<string> checked_line;
};

struct kmer_smoothed {
    vector<float> signal_smoothed;
    vector<float> distance_vector;
    string id_kmer;
};

struct signal_t{
    string ID;
    vector<vector<float>> signals;

    template <class B>
    void serialize(B& buf) const {
        buf << ID << signals;
    }

    template <class B>
    void parse(B& buf) {
        buf >> ID  >> signals;
    }                      
};
