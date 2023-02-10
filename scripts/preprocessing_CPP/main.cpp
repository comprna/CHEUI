/**
 * @author Wenjing Xue
 * @brief
 * @version 0.1
 * @date 2022-11-10
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <future>
#include <mutex>
#include <set>
#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>
#include <bits/stdc++.h>

#include "argagg.hpp"
#include "kmer.hpp"
#include "PicklingTools/chooseser.h"

using namespace std;
const int RUN = 32; // Initialising the RUN to get chunks

/**
 * @brief
 * A simple function to split a giving string with a delimiter
 */
vector<string> split_string(string str, char delimiter)
{
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while (getline(ss, tok, delimiter))
    {
        internal.push_back(tok);
    }

    return internal;
}

/**
 * @brief
 * A function to read the csv file
 */
vector<kmer_models> read_csv(string file, char delimiter)
{
    vector<kmer_models> kmers;
    kmer_models kmer;
    ifstream infile(file);
    string line;
    string word;
    int dataID = 1;
    int num = 0;

    getline(infile, line);

    while (getline(infile, line))
    {
        stringstream str(line);

        while (getline(str, word, delimiter))
        {
            if (dataID == 1)
            {
                kmer.kmers = word;
                ++dataID;
            }
            else if (dataID == 2)
            {
                kmer.mean = stof(word);
                ++dataID;
            }
            else
            {
                kmer.stdv = stof(word);
                dataID = 1;
                kmers.push_back(kmer);
            }
        }
    }

    return kmers;
}

/**
 * @brief
 * split the large nanopolish file into several small files
 * to make it run on multiple threads
 */
set<string> split_file(string nanopolish_path, int num_file, string dir_out, char symbol)
{
    set<string> filepath;
    ifstream infile(nanopolish_path);
    string line;
    getline(infile, line);

    // check the header
    string headerline = line;
    vector<string> headers = split_string(headerline, '\t');

    if (headers[headers.size() - 1] != "samples")
    {
        throw "Error: nanopolish samples are not found, please run nanopolish with flag --samples";
    }

    // get the column index
    int line_count = 0;
    while (getline(infile, line))
    {
        if (line != "\n")
        {
            ++line_count;
        }
    }
    ++line_count;

    infile.clear();
    infile.seekg(0, ios::beg);
    int counter_reads = 0;
    int chunk_size = line_count / num_file;
    int i = 0;
    bool header = false;

    vector<string> N_lines_print;
    getline(infile, line);
    ofstream f;

    while (getline(infile, line))
    {
        ++counter_reads;
        if (counter_reads < chunk_size)
        {
            if (!header)
            {
                N_lines_print.push_back(headerline);
                header = true;
            }

            N_lines_print.push_back(line);

            if (N_lines_print.size() > 200000)
            {
                string file = dir_out + "/temp_" + to_string(i) + ".tmp";
                filepath.insert(file);
                f.open(file, ios::app);

                for (string s : N_lines_print)
                {
                    f << s << "\n";
                }
                f.close();
                N_lines_print.clear();
            }
        }
        else
        {
            if (N_lines_print.size() > 0)
            {
                string file = dir_out + "/temp_" + to_string(i) + ".tmp";
                filepath.insert(file);
                f.open(file, ios::app);

                for (string s : N_lines_print)
                {
                    f << s << "\n";
                }
                f.close();
                N_lines_print.clear();
            }

            if (split_string(line, '\t')[2].find(symbol) != string::npos)
            {
                string file = dir_out + "/temp_" + to_string(i) + ".tmp";
                filepath.insert(file);
                f.open(file, ios::app);
                f << line << "\n";
                f.close();
            }
            else
            {
                string file = dir_out + "/temp_" + to_string(i) + ".tmp";
                filepath.insert(file);
                f.open(file, ios::app);
                f << line << "\n";
                f.close();
                ++i;
                counter_reads = 0;
                header = false;
            }
        }
    }

    string file = dir_out + "/temp_" + to_string(i) + ".tmp";
    filepath.insert(file);
    f.open(file, ios::app);

    for (string s : N_lines_print)
    {
        f << s << "\n";
    }
    f.close();
    N_lines_print.clear();

    return filepath;
}

int get_index(vector<string> v, string s)
{
    int index = -1;
    auto it = find(v.begin(), v.end(), s);
    if (it != v.cend())
    {
        index = distance(v.begin(), it);
    }

    return index;
}

vector<pair<int, kmers_parse>> _recycle_kmers(vector<pair<int, kmers_parse>> &kmer_dict, char symbol)
{
    vector<pair<int, kmers_parse>> keep_kmers_dict;
    int size = kmer_dict.size();
    vector<int> reverse_list;

    for (auto k : kmer_dict)
    {
        reverse_list.emplace(reverse_list.begin(), k.first);
    }

    for (int i = 0; i < size - 1; ++i)
    {
        if ((reverse_list[i] - 1) == reverse_list[i + 1])
        {
            keep_kmers_dict.push_back(kmer_dict[size - i - 1]);
        }
        else
        {
            break;
        }
    }

    int i = 0;
    while (i < keep_kmers_dict.size())
    {
        if (keep_kmers_dict[i].second.kmers.back() != symbol)
        {
            keep_kmers_dict.erase(keep_kmers_dict.begin() + i);
        }
        else
        {
            ++i;
        }
    }
    return keep_kmers_dict;
}

/**
 * @brief
 * Do a regex search against all defined regexes and
 * return the key and match result of the first matching regex
 */
vector<string> _check_line(string line, int contig_idx, int position_idx, int reference_kmer_idx, int model_kmer_idx, int samples_idx, char symbol)
{
    vector<string> line_split;
    line_split = split_string(line, '\t');

    // check model and reference are the same
    if (line_split[reference_kmer_idx] != line_split[model_kmer_idx])
    {
        return {};
    }

    // check the model is not NNNNN
    if (line_split[model_kmer_idx] == "NNNNN")
    {
        return {};
    }

    // check if there is any A in the model
    if (line_split[model_kmer_idx].find(symbol) == string::npos)
    {
        return {};
    }

    return line_split;
}

// This function sorts array from left index to
// to right index which is of size atmost RUN
template <typename T>
void insertionSort(vector<T> &arr, int left, int right)
{
    for (int i = left + 1; i <= right; ++i)
    {
        T temp = arr[i];
        int j = i - 1;
        while (j >= left && arr[j] > temp)
        {
            arr[j + 1] = arr[j];
            --j;
        }
        arr[j + 1] = temp;
    }
}

template <typename It, typename T>
void insertion_sort(vector<pair<It, T>> &arr, int left, int right)
{
    for (int i = left + 1; i <= right; ++i)
    {
        pair<It, T> temp = arr[i];
        int j = i - 1;
        while (j >= left && arr[j].first > temp.first)
        {
            arr[j + 1] = arr[j];
            --j;
        }
        arr[j + 1] = temp;
    }
}

// Merge function merges the sorted runs
template <typename T>
void merge(vector<T> &arr, int l, int m, int r)
{

    // Original array is broken in two parts left and right array
    int len1 = m - l + 1, len2 = r - m;
    T left[len1], right[len2];
    for (int i = 0; i < len1; ++i)
        left[i] = arr[l + i];
    for (int i = 0; i < len2; ++i)
        right[i] = arr[m + 1 + i];

    int i = 0;
    int j = 0;
    int k = l;

    // After comparing, we merge those two array in larger sub array
    while (i < len1 && j < len2)
    {
        if (left[i] <= right[j])
        {
            arr[k] = left[i];
            ++i;
        }
        else
        {
            arr[k] = right[j];
            ++j;
        }
        ++k;
    }

    // Copy remaining elements of left, if any
    while (i < len1)
    {
        arr[k] = left[i];
        ++k;
        ++i;
    }

    // Copy remaining element of right, if any
    while (j < len2)
    {
        arr[k] = right[j];
        ++k;
        ++j;
    }
}

template <typename It, typename T>
void _merge(vector<pair<It, T>> &arr, int l, int m, int r)
{

    // Original array is broken in two parts left and right array
    int len1 = m - l + 1, len2 = r - m;
    pair<It, T> left[len1], right[len2];
    for (int i = 0; i < len1; ++i)
        left[i] = arr[l + i];
    for (int i = 0; i < len2; ++i)
        right[i] = arr[m + 1 + i];

    int i = 0;
    int j = 0;
    int k = l;

    // After comparing, we merge those two array in larger sub array
    while (i < len1 && j < len2)
    {
        if (left[i].first <= right[j].first)
        {
            arr[k] = left[i];
            ++i;
        }
        else
        {
            arr[k] = right[j];
            ++j;
        }
        ++k;
    }

    // Copy remaining elements of left, if any
    while (i < len1)
    {
        arr[k] = left[i];
        ++k;
        ++i;
    }

    // Copy remaining element of right, if any
    while (j < len2)
    {
        arr[k] = right[j];
        ++k;
        ++j;
    }
}

template <typename It, typename T>
void tim_sort(vector<pair<It, T>> &arr)
{
    int n = arr.size();

    // Sort individual subarrays of size RUN
    for (int i = 0; i < n; i += RUN)
        insertion_sort(arr, i, min((i + RUN - 1), (n - 1)));

    // Start merging from size RUN (or 32).
    for (int size = RUN; size < n; size += size)
    {

        // pick starting point of left sub array.
        // We are going to merge
        // arr[left..left+size-1]
        // and arr[left+size, left+2*size-1]
        // After every merge, we
        // increase left by 2*size
        for (int left = 0; left < n; left += (size + size))
        {

            // find ending point of left sub array mid+1 is starting point of right sub array
            int mid = left + size - 1;
            int right = min((left + size + size - 1), (n - 1));

            // merge sub array arr[left.....mid] &
            // arr[mid+1....right]
            if (mid < right)
                _merge(arr, left, mid, right);
        }
    }
}

// Iterative Timsort function to sort the
// array[0...n-1] (similar to merge sort)
// template <class It> void timSort(vector<It> arr, int n)
template <typename T>
void timSort(vector<T> &arr)
{
    int n = arr.size();
    // Sort individual subarrays of size RUN
    for (int i = 0; i < n; i += RUN)
        insertionSort(arr, i, min((i + RUN - 1), (n - 1)));

    // Start merging from size RUN (or 32).
    // It will merge
    // to form size 64, then 128, 256
    // and so on ....
    for (int size = RUN; size < n; size += size)
    {

        // pick starting point of
        // left sub array. We
        // are going to merge
        // arr[left..left+size-1]
        // and arr[left+size, left+2*size-1]
        // After every merge, we
        // increase left by 2*size
        for (int left = 0; left < n; left += (size + size))
        {

            // find ending point of
            // left sub array
            // mid+1 is starting point
            // of right sub array
            int mid = left + size - 1;
            int right = min((left + size + size - 1), (n - 1));

            // merge sub array arr[left.....mid] &
            // arr[mid+1....right]
            if (mid < right)
                merge(arr, left, mid, right);
        }
    }
}

/**
 * @brief
 *
 */
kmers_parse_r _parse_kmers(vector<string> &checkedLine, int contig_idx, int position_idx, int reference_kmer_idx, int model_kmer_idx, int samples_idx, ifstream &infile, vector<pair<int, kmers_parse>> &parsed_kmer, int &counter, char symbol)
{
    vector<pair<int, kmers_parse>> kmer_lines;
    string line;
    int pos1;

    if (!parsed_kmer.empty())
    {
        tim_sort(parsed_kmer);
        kmer_lines = parsed_kmer;
        pos1 = kmer_lines[0].first;
    }
    else
    {
        // kmer_lines.clear();
        pos1 = stoi(checkedLine[position_idx]);
    }

    vector<int> positions = {pos1, pos1 + 1, pos1 + 2, pos1 + 3, pos1 + 4};

    while (kmer_lines.size() < 5)
    {
        if (!checkedLine.empty())
        {
            vector<float> samples;
            for (string i : split_string(checkedLine[samples_idx], ','))
            {
                samples.push_back(stof(i));
            }

            // find the val in kmer_lines.keys(), check whether it is in the dict
            int val = stoi(checkedLine[position_idx]);
            auto it = find_if(kmer_lines.begin(), kmer_lines.end(), [&](const pair<int, kmers_parse> &element)
                              { return element.first == val; });
            if (it != kmer_lines.end())
            {
                int index = distance(kmer_lines.begin(), it);

                kmer_lines[index].second.samples.insert(kmer_lines[index].second.samples.end(), samples.begin(), samples.end());
            }
            else
            {
                pair<int, kmers_parse> k;
                k.first = stoi(checkedLine[position_idx]);
                k.second = kmers_parse{checkedLine[model_kmer_idx], samples, checkedLine[contig_idx], checkedLine[3]};
                kmer_lines.push_back(k);
            }
        }

        if (getline(infile, line))
        {
            ++counter;
        }
        else
        {
            return kmers_parse_r{{}, counter, {}};
        }

        checkedLine = _check_line(line, contig_idx, position_idx, reference_kmer_idx, model_kmer_idx, samples_idx, symbol);

        if (!checkedLine.empty() && ((stoi(checkedLine[1]) > (pos1 + 4)) || (stoi(checkedLine[1]) < pos1)))
        {
            break;
        }
        else
        {
            continue;
        }
    }

    // check if it forms a 9mer
    if (kmer_lines.size() < 5)
    {
        return kmers_parse_r{{}, counter, checkedLine};
    }

    tim_sort(kmer_lines);
    vector<int> keys;
    for (auto k : kmer_lines)
    {
        keys.push_back(k.first);
    }

    if (keys == positions)
    {
        return kmers_parse_r{kmer_lines, counter, checkedLine};
    }
    else
    {
        // grab the index of the mistmatch between the expected kmers and the
        // ones extracted without fail
        vector<int> matches;
        for (int i = 0; i < keys.size(); ++i)
        {
            if (keys[i] == positions[i])
            {
                matches.push_back(i);
            }
        }

        for (int i = matches.size() - 1; i >= 0; --i)
        {
            kmer_lines.erase(kmer_lines.begin() + matches[i]);
        }

        // check if numbers are consecutive and take the last consecutive ones ending in A/C
        kmer_lines = _recycle_kmers(kmer_lines, symbol);

        return kmers_parse_r{kmer_lines, counter, checkedLine};
    }
}

float median(vector<float> array)
{
    if (array.empty())
    {
        return 0;
    }
    else
    {
        int size = array.size();
        float median;
        vector<float> temp = array;
        // stable_sort(temp.begin(), temp.end());
        // sort(temp.begin(), temp.end());
        timSort(temp);
        if (size % 2 == 0)
        {
            return (temp[size / 2 - 1] + temp[size / 2]) / 2.0;
        }
        else
        {
            return temp[size / 2];
        }
    }
}

/**
 * @brief
 * This function top an array until some specific lenght
 */
vector<float> top_median(vector<float> array)
{
    int size = array.size();
    float m = median(array);
    for (int i = size; i < 20; ++i)
    {
        array.push_back(m);
    }

    return array;
}

vector<float> slicing(vector<float> &arr, int x, int y)
{
    // Starting and Ending iterators
    auto start = arr.begin() + x;
    auto end = arr.begin() + y + 1;

    // To store the sliced vector
    vector<float> result(y - x + 1);

    // Copy vector using copy function()
    copy(start, end, result.begin());

    // Return the final sliced vector
    return result;
}

/**
 * @brief
 * smmoth the signal
 */
vector<float> smooth_event(vector<float> raw_signal)
{
    vector<float> raw_signal_events;

    if (raw_signal.size() < 20)
    {
        vector<float> event = top_median(raw_signal);
        for (float val : event)
        {
            // val = (int)(val * 1000.0) / 1000.0;
            raw_signal_events.push_back(val);
        }
    }
    else
    {
        int division = raw_signal.size() / 20;
        vector<float> new_event;
        for (int i = 0; i < raw_signal.size(); i += division)
        {
            float m = median(slicing(raw_signal, i, i + division));
            new_event.push_back(m);

            if (new_event.size() == 20)
            {
                break;
            }
        }

        if (new_event.size() < 20)
        {
            new_event = top_median(new_event);
        }

        for (float val : new_event)
        {
            // val = (int)(val * 1000.0) / 1000.0;
            raw_signal_events.push_back(val);
        }
    }

    return raw_signal_events;
}

vector<float> make_expected(string kmer, vector<kmer_models> &model_kmer)
{
    vector<float> expected_signal;
    for (int i = 0; i < 5; ++i)
    {
        string kmer_5 = kmer.substr(i, 5);
        for (auto model : model_kmer)
        {
            if (model.kmers == kmer_5)
            {
                vector<float> signal_add(20, model.mean);
                expected_signal.insert(expected_signal.end(), signal_add.begin(), signal_add.end());
            }
        }
    }
    return expected_signal;
}

vector<float> distance_calculator(vector<float> signal_expected, vector<float> event_smoothed)
{
    vector<float> vector_distance;
    for (int i = 0; i < signal_expected.size(); ++i)
    {
        float val = abs(signal_expected[i] - event_smoothed[i]);
        // val = (int)(val * 1000.0) / 1000.0;
        vector_distance.push_back(val);
    }

    return vector_distance;
}

/**
 * @brief
 * Smooth the signals to fix the length
 */
kmer_smoothed _smooth_kmer(vector<pair<int, kmers_parse>> &parsed_kmer, vector<kmer_models> &model_kmer)
{
    vector<pair<int, kmers_parse>> sorted = parsed_kmer;

    string kmer_name = sorted[0].second.kmers + sorted[1].second.kmers.back() + sorted[2].second.kmers.back() + sorted[3].second.kmers.back() + sorted[4].second.kmers.back();

    string id_kmer = parsed_kmer[0].second.contig + "_" + to_string(sorted[0].first) + "_" + kmer_name + "_" + sorted[0].second.name;

    vector<float> signal_smoothed;
    for (auto p : sorted)
    {
        vector<float> event = p.second.samples;
        // smooth the event
        vector<float> event_smoothed = smooth_event(event);
        // add the event to the signal
        signal_smoothed.insert(signal_smoothed.end(), event_smoothed.begin(), event_smoothed.end());
    }

    // create an expected signal according to the kmer
    vector<float> expected_smoothed = make_expected(kmer_name, model_kmer);

    // claculate distance between expected and actual signal
    vector<float> distance_vector = distance_calculator(expected_smoothed, signal_smoothed);

    return kmer_smoothed{signal_smoothed, distance_vector, id_kmer};
}

/**
 * @brief
 * combine signals and distance vectors
 */
Arr _combine_vectors(vector<float> smooth_signal, vector<float> smooth_distance)
{
    Arr combined(100);
    for (int i = 0; i < smooth_signal.size(); ++i)
    {
        Arr c(2);
        c.append(smooth_signal[i]);
        c.append(smooth_distance[i]);
        combined.append(c);
    }

    return combined;
}

/**
 * @brief
 *
 * @param file
 */
void parse_nanopolish(string file, char symbol, vector<kmer_models> &model_kmer, string name_out)
{
    // vector<kmer_models> model_kmer = read_csv(model_kmer_path, ',');
    int counter = 0;
    vector<pair<int, kmers_parse>> parsed_kmer;
    vector<string> stored_line;
    vector<int> signals_IDs;
    vector<string> checked_line;

    ifstream infile(file);
    string line;
    getline(infile, line);

    // check the header
    string headerline = line;
    vector<string> header = split_string(headerline, '\t');

    if (header[header.size() - 1] != "samples")
    {
        throw "Error: nanopolish samples are not found, please run nanopolish with flag --samples";
    }

    // get the column index
    int contig_idx = get_index(header, "contig");
    int position_idx = get_index(header, "position");
    int reference_kmer_idx = get_index(header, "reference_kmer");
    int model_kmer_idx = get_index(header, "model_kmer");
    int samples_idx = get_index(header, "samples");

    if (contig_idx == -1 || position_idx == -1 || reference_kmer_idx == -1 || model_kmer_idx == -1 || samples_idx == -1)
    {
        throw "Some nanopolish columns are not found\n";
    }

    // The EOF char is an empty string
    while (!infile.eof())
    {
        if (stored_line.empty())
        {
            if (getline(infile, line).fail())
            {
                break;
            }
            ++counter;

            // check line is fordward and does not have NNNNN in the model
            checked_line = _check_line(line, contig_idx, position_idx, reference_kmer_idx, model_kmer_idx, samples_idx, symbol);
        }
        else
        {
            checked_line = stored_line;
            stored_line.clear();
        }

        if (!checked_line.empty())
        {
            if ((checked_line[model_kmer_idx].back() == symbol) || (!parsed_kmer.empty()))
            {
                kmers_parse_r kmers_p = _parse_kmers(checked_line, contig_idx, position_idx, reference_kmer_idx, model_kmer_idx, samples_idx, infile, parsed_kmer, counter, symbol);
                parsed_kmer = kmers_p.kmers_lines;
                stored_line = kmers_p.checked_line;

                //  if parsed kmer fail the below code does not get executed and new_parser_kmers never get recycle
                if (!parsed_kmer.empty())
                {
                    if (parsed_kmer.size() < 5)
                    {
                        continue;
                    }

                    kmer_smoothed kmer_s = _smooth_kmer(parsed_kmer, model_kmer);
                    vector<float> smooth_signal = kmer_s.signal_smoothed;
                    vector<float> smooth_distance = kmer_s.distance_vector;
                    string ID = kmer_s.id_kmer;

                    Arr combined_signals = _combine_vectors(smooth_signal, smooth_distance);

                    Tab signals_IDs;
                    signals_IDs[ID] = combined_signals;

                    DumpValToFile(signals_IDs, name_out, SERIALIZE_P2);
                    signals_IDs = {};

                    vector<string> IDs = split_string(ID, '_');
                    string s = IDs[IDs.size() - 2].substr(5);
                    int index_ = s.find_first_of(symbol);

                    if (index_ != -1)
                    {
                        index_ = index_ + 1;

                        vector<pair<int, int>> sorted;
                        vector<int> delete_positions;
                        for (int i = 0; i < parsed_kmer.size(); ++i)
                        {
                            sorted.push_back(pair<int, int>{parsed_kmer[i].first, i});
                        }
                        tim_sort(sorted);

                        for (int i = 0; i < index_; ++i)
                        {
                            delete_positions.push_back(sorted[i].second);
                        }

                        timSort(delete_positions);

                        for (int i = delete_positions.size() - 1; i >= 0; --i)
                        {
                            parsed_kmer.erase(parsed_kmer.begin() + delete_positions[i]);
                        }
                                              
                    }
                    else
                    {
                        parsed_kmer.clear();
                        // recover the information with the kmers
                        continue;
                    }
                }
            }
        }

        if (counter % 500000 == 0)
        {
            cout << counter << " processd lines" << endl;
        }
    }
}

/**
 * @brief
 * Main function
 */
int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Run with mode: ./CHEUI preprocess <m6A/m5C>" << std::endl;
        return EXIT_FAILURE;
    }

    argagg::parser argparser{{
        {"help", {"-h", "--help"}, "shows this help message", 0},
        {"input", {"-i", "--input-nanopolish"}, "Nanopolish file. <out.txt>", 1},
        {"model", {"-m", "--kmer-model"}, "file containing all the expected signal k-mer means", 1},
        {"output", {"-o", "--out-dir"}, "output folder (default: .)", 1},
        {"threads", {"-n", "--CPU"}, "Number of CPUs (threads) to use (default: 1)", 1},
        {"suffix", {"-s", "--suffix-name"}, "name to use for output files", 1},
        {"type1", {"--m6A"}, "select the preprocess type", 0},
        {"type2", {"--m5C"}, "select the preprocess type", 0},
    }};

    argagg::parser_results args;
    try
    {
        args = argparser.parse(argc, argv);
    }
    catch (const exception &e)
    {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }

    if (args["help"])
    {
        cerr << argparser;
        return EXIT_SUCCESS;
    }

    if (!args["input"])
    {
        cerr << "ERROR: No input file provided" << endl;
        cerr << argparser;
        return EXIT_FAILURE;
    }

    char symbol;
    bool type1 = args["type1"];
    bool type2 = args["type2"];
    if (type1 && !type2)
    {
        symbol = 'A';
    }
    else if (type2 && !type1)
    {
        symbol = 'C';
    }
    else
    {
        cerr << "Error: Please select a correct data type (m6A/m5C)" << endl;
        return EXIT_FAILURE;
    }

    int n_threads = args["threads"].as<int>(1);

    string model_kmer_path = args["model"].as<string>("");

    string outdir = args["output"].as<string>(".");
    string name_out;
    if (args["suffix"])
    {
        name_out = outdir + "/" + args["suffix"].as<string>("") + "_signals+IDS.p";
    }
    else
    {
        vector<string> path = split_string(args["input"].as<string>(""), '/');
        string p = path[path.size() - 1];
        p = p.substr(0, p.size() - 4);
        name_out = outdir + "/" + p + "_signals+IDS.p";
    }

    if (mkdir(outdir.c_str(), 0777) == -1)
    {
        if (remove(name_out.c_str()) == 0)
        {
            cerr << "\nWARNING: Signal file already exist, deleting the previous generated signals+IDs files \n"
                 << endl;
        }
    }

    if (remove(name_out.c_str()) == 0)
    {
        cerr << "\nWARNING: Signal file already exist, deleting the previous generated signals+IDs files \n"
             << endl;
    }

    // if threads > 1, split the input file and process the tmp file in parallel method
    // otherwise, direct process the original input file
    if (n_threads > 1)
    {
        set<string> filepath;
        try
        {
            filepath = split_file(args["input"].as<string>(""), n_threads, args["output"].as<string>(""), symbol);
        }
        catch (const char *c)
        {
            cerr << c;
            return EXIT_FAILURE;
        }

        vector<future<int>> tasks;
        vector<kmer_models> model_kmer;
        for (string file : filepath)
        {
            tasks.emplace_back(std::async(std::launch::async, [file, &symbol, &model_kmer_path, &name_out]
                                          {
                vector<kmer_models> model_kmer = read_csv(model_kmer_path, ',');
                try {
                    parse_nanopolish(file, symbol, model_kmer, name_out);
                    return EXIT_FAILURE;
                }
                catch (const char* c) {
                    cerr << c;
                    return EXIT_FAILURE;
                } }));
        }

        for (auto &&task : tasks)
        {
            task.get();
        }

        for (string file : filepath)
        {
            remove(file.c_str());
        }

        // #ifndef _WIN32
        //         string command = "rm ";
        //         string dir = outdir + "/*.tmp";
        //         system(command.append(dir).c_str());
        // #endif
    }
    else
    {
        vector<kmer_models> model_kmer = read_csv(model_kmer_path, ',');
        try
        {
            parse_nanopolish(args["input"].as<string>(""), symbol, model_kmer, name_out);
        }
        catch (const char *c)
        {
            cerr << c;
            return EXIT_FAILURE;
        }
    }

    return 0;
}
