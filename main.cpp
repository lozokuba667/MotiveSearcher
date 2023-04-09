#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <utility>
#include <vector>
#include <tuple>


using namespace std;

class Vert {
public:
     string substring;
     int sequence_number;
     int position;
};




tuple<vector<vector<char>>, vector<vector<int>>> readFiles() {
    /*
 * PARAMETRY: brak
 * ZWRACA: krotkę wektorów z wczytanymi instancjami plików sekwencji i jakości
 * Funkcja służąca do wczytania testowych instancji
 * Iteruje plik linia po linii, pomija identyfikatory sekwencji, sprowadza linie do postaci strumienia
 * W zależności od obecnego pliku, wkłada dane do odpowiedniego wektora
 */
    vector<vector<char>> seq(5);
    vector<vector<int>> quality(5);
    vector<string> file_names = {"sequence.fasta", "sample.qual"};

    for (const auto &name: file_names) {
        fstream file(name);
        int line_counter = 0;
        string sinleline;

        while (getline(file, sinleline)) {
            if (sinleline.find('>') != string::npos) { continue; }
            istringstream streamline(sinleline);
            if (name == "sequence.fasta") {
                seq[line_counter] = vector<char>(istream_iterator<char>(streamline), istream_iterator<char>());
            } else {
                quality[line_counter] = vector<int>(istream_iterator<int>(streamline), istream_iterator<int>());
            }
            line_counter++;
        }
        file.close();
    }
    return {seq, quality};
}


void showContent(vector<vector<char>> sequences, vector<vector<int>> qualities) {
    cout << "\n";
    for (int i = 0; i < sequences.size(); i++) {
        cout << "Seq " << i + 1 << ": ";
        for (char singleSeq: sequences[i]) {
            cout << singleSeq;
        }
        cout << "\n";
        cout << "Quality measure for Seq " << i + 1 << ": " << endl;
        for (int singleQuality: qualities[i]) {
            cout << singleQuality << " ";
        }
        cout << "\n";
        cout << "\n";
    }
}
vector<vector<int>> rememeberPositions(vector<vector<int>> quality, int limit){
    vector<vector<int>> positions(5);
    for(int i=0; i<quality.size(); i++){
        for(int j=0; j<quality[i].size();j++){
            if(quality[i][j] >=limit){positions[i].push_back(j + 1);}
        }
    }
    return positions;
}

void eraseBelowQuality(vector<vector<char>> &seq, vector<vector<int>> &qual, int limit) {
    for (int i = 0; i < qual.size(); i++) {
        for (int j = 0; j < qual[i].size(); j++) {
            if (qual[i][j] < limit) {
                qual[i].erase(qual[i].begin() + j);
                seq[i].erase(seq[i].begin() + j);
                j=j-1;
            }
        }
    }
}

vector<vector<int>> qualityBasedTrimmingWrapper(vector<vector<char>> &seq, vector<vector<int>> &qual) {
    int limit;
    cout << "Enter limit for quality: ";
    cin >> limit;
    vector<vector<int>> onlyGoodPositions = rememeberPositions(qual,limit);
    eraseBelowQuality(seq,qual,limit);
    return onlyGoodPositions;
}

vector<Vert> createListOfVerticles(vector<vector<char>> seq, vector<vector<int>> pos, int substring_size){
    vector<Vert> verts;
    for(int i=0; i<seq.size(); i++) {
        for (int j = 0; j < seq[i].size() - substring_size + 1; j++) {
            Vert temp;
            for (int z = j; z < substring_size + j; z++) {
                temp.substring.push_back(seq[i][z]);
            }
            temp.position = pos[i][j];
            temp.sequence_number = i;
            verts.push_back(temp);
        }
    }

    return verts;
}




bool isProperForConnect(Vert v1, Vert v2, int substring_size){
    if(v1.substring == v2.substring){
        if(v1.sequence_number != v2.sequence_number){
            int higher_pos = max(v1.position, v2.position);
            int lower_pos = min(v1.position, v2.position);
            if(higher_pos - lower_pos <= 10*substring_size){
                return true;
            }
        }
    }
    return false;
}


vector<pair<Vert,vector<Vert>>> createGraph(vector<vector<char>> sequences, vector<vector<int>> positons){
    int substring_size;
    cout << "Enter desired substring size: ";
    cin >> substring_size;

    vector<Vert> verticles = createListOfVerticles(std::move(sequences), std::move(positons),substring_size);
    vector<pair<Vert,vector<Vert>>> graph;
    vector<Vert> rememberedConnection;
    for(int i=0; i<verticles.size()-1; i++){
        for(int j=i; j<verticles.size(); j++){
            if(isProperForConnect(verticles[i], verticles[j], substring_size)){
                rememberedConnection.push_back(verticles[j]);
            }
        }
        graph.emplace_back(verticles[i], rememberedConnection);
        rememberedConnection.clear();
    }
    return graph;
}

void showGraph(vector<pair<Vert,vector<Vert>>> graph){
    cout << "\n Pseudo-Successor List\n";
    for(auto item:graph){
        cout << "[" << item.first.substring << ", S:" << item.first.sequence_number << ", P:" << item.first.position << "]" << " ---> ";
        for(auto subitem: item.second){
            cout << "(S:" << subitem.sequence_number << ", P:" << subitem.position << "),";
        }
        cout<<"\n";
    }
}

bool isSingleSequenceConnected(vector<Vert> current){
    vector<int> seq_numbers;
    for(auto item:current){
        seq_numbers.push_back(item.sequence_number);
    }
    if(adjacent_find(seq_numbers.begin(),seq_numbers.end()) == seq_numbers.end()){return true;}
    return false;
}

int reassureConnections(vector<Vert> verts, int index, vector<pair<Vert,vector<Vert>>> graph){
    int points = 0;
    vector<Vert> sliced = verts;
    for(auto successor:verts){
        if(successor.sequence_number == 4){
            points+=1;
            return points;
        } else {
            sliced = {sliced.begin() + 1, sliced.end()};
            for(int i=index; i<graph.size();i++){
                if(successor.sequence_number == graph[i].first.sequence_number && successor.position == graph[i].first.position){
                    for(int j=0; j<graph[i].second.size(); j++){
                        if(sliced[j].sequence_number == graph[i].second[j].sequence_number && sliced[j].position == graph[i].second[j].position) {points+=1;}
                        else {points-=1;}
                    }
                }
            }
        }
    }
    return points;
}

void showMotive(pair<Vert,vector<Vert>> motive){
    cout << "\nMOTIVE: " << motive.first.substring << "\n";
    cout << "Found in sequence {" << motive.first.sequence_number << "} on position {" << motive.first.position << "}\n";
    for(auto item:motive.second){
        cout << "Found in sequence {" << item.sequence_number << "} on position {" << item.position << "}\n";

    }
}

void findMotive(vector<pair<Vert,vector<Vert>>> graph){
    int connection_points = 0;
    pair<Vert,vector<Vert>> motive;
    for(int i=0; i<graph.size(); i++){
        if(graph[i].second.size() == 4 && graph[i].first.sequence_number == 0){
            if(isSingleSequenceConnected(graph[i].second)){
                connection_points = reassureConnections(graph[i].second, i+1, graph);
                if(connection_points == 7) {motive = graph[i]; break;}
                else{connection_points = 0;}
            }
        } else if(graph[i].first.sequence_number != 0){
            cout << "Clear Motive is not present in provided data\n";
            exit(0);
        }
    }
    showMotive(motive);
}

int main() {
    auto[instance_seq, instance_qual] = readFiles();
    cout << "\nREAD INSTANCE (before trimming):\n";
    showContent(instance_seq, instance_qual);
    vector<vector<int>> properPositions = qualityBasedTrimmingWrapper(instance_seq,instance_qual);
    cout << "\nREAD INSTANCE (after trimming):\n";
    showContent(instance_seq,instance_qual);
    vector<pair<Vert,vector<Vert>>> graph = createGraph(instance_seq,properPositions);
    showGraph(graph);
    findMotive(graph);
    return 0;
}




