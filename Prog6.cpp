#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <queue>

using namespace std;

#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

const float earthradius = 3963.0;     // [miles]
const float distance_accuracy = 5.0;  // [miles]

const int national_minpop = 1000000;

const float national_dist = 150.0;    // [miles]
const float regional_dist = 100.0;    // [miles]

const float local_maxdist = 50.0;     // [miles]

const float plane_speed = 520.0;      // [mph]
const float truck_speed = 65.0;       // [mph]

enum city_t { LOCAL, METRO, REGIONAL, NATIONAL };
enum color_t { WHITE, BLACK };

// City labels
// NATIONAL: national_minpop < metro_pop
//           maxpop < national_dist
// REGIONAL: maxpop < regional_dist
// LOCAL:    otherwise

// Edge connectivity (bidirectional)
// NATIONAL: all NATIONAL cities 
// REGIONAL: three nearest NON-LOCAL cities
// LOCAL:    five nearest NON-LOCAL cities
//           all LOCAL cities < local_maxdist

class city { 
    friend istream& operator>>(istream &, city &);
    friend ostream& operator<<(ostream &, const city &);
  
  public:
    
    /*city(const string &name, float latitude, float longitude, int population, city_t type) {
      this->name = name;
      this->latitude = latitude;
      this->longitude = longitude;
      this->population = population;
      this->type = type;
    }*/

    string getName() const { return name; }
    float getLatitude() const { return latitude; }
    float getLongitude() const { return longitude; }
    int getPopulation() const { return population; }
    city_t getType() const { 
      return type;
    }

    void setName(const std::string& name) { this->name = name; }
    void setLatitude(float latitude) { this->latitude = latitude; }
    void setLongitude(float longitude) { this->longitude = longitude; }
    void setPopulation(int population) { this->population = population; }
    void setType(city_t type) { this->type = type; }

    int internal_index;
    static int w;
    
  private:
    string name;
    float latitude;
    float longitude;
    int population;
    city_t type;
};

int city::w = 0;

class matrix {

  private:
    vector<vector<float>> data;
    int size;

  public:
    matrix (int n) : size(n) {
      data.resize(n);
      for (int i = 0; i < n; i++) {
        data[i].resize(i + 1, 0.0);
      }
    }

    float getValue(int row, int col) const {
      //ensure row is greater than or equal to col
      if (row < col) {
        std::swap(row, col);
      }
      return data[row][col];
    }

    void setValue(int row, int col, float value) {
      if (row < col) {
        std::swap(row, col);
      }
      data[row][col] = value;
    }

    int getSize() const {
      return size;
    }
 };

istream & operator>>(istream &in, city &c) {

    //get data from each line
    string line;
    if (getline(in, line)) {

      stringstream ss(line);
      string city, state, type_str;
      float latitude, longitude;
      int population;

      getline(ss, city, ',');
      getline(ss, state, ',');
      getline(ss, type_str, ',');
      ss >> latitude;
      ss.ignore();
      ss >> longitude;
      ss.ignore();
      ss >> population;

      //append state to city name
      city += "_" + state;

      
      if (type_str == "METRO") {
        c.setType(METRO);
      }
      else {
        c.setType(LOCAL);
      }

      c.name = city;
      c.latitude = latitude;
      c.longitude = longitude;
      c.population = population;

      //set city w for max width
      if((int)c.name.length() >= city::w) {
        city::w = c.name.length();
      }
    }
    return in;
}
ostream & operator<<(ostream &out, const city &c) {

    int dots = city::w - c.name.length() + 3;
    string type_str;
    //change enum to the string name for types
    switch (c.getType())
    {
    case 0:
      type_str = "LOCAL";
      break;
    case 1:
      type_str = "METRO";
      break;
    case 2:
      type_str = "REGIONAL";
      break;
    case 3:
      type_str = "NATIONAL";
      break;
    }

    out << c.name;
    for (int i = 0; i < dots; i++) {
      out << ".";
    }

    out << "  " << setw(8) << left << type_str << "  " << setw(8) << right << c.population << "  " << setw(7) << right << fixed << setprecision(2)
        << c.latitude << "  " << setw(7) << fixed << setprecision(2) << c.longitude;

    return out;
}


bool compare_population(const city& a, const city& b) {
  return a.getPopulation() > b.getPopulation();
}

void create_vertex_table(istream &in, vector<city> &vertex_table) { 

    city::w = 0;
    city c;

    //push back cities
    while (in >> c) {
      vertex_table.push_back(c);
    }

    //sort by highest population
    std::sort(vertex_table.begin(), vertex_table.end(), compare_population);

    for (int i=0; i < vertex_table.size(); i++) {
      vertex_table[i].internal_index = i;
    }

}

void update_vertex_table(vector<city> &vertex_table, matrix &dist_table) {

    //Rule 1: relabel metro cities based on population
    for (city& city : vertex_table) {
      if (city.getType() == METRO) {
        if (city.getPopulation() >= 1000000) {
          city.setType(NATIONAL);
        } else {
          city.setType(REGIONAL);
        }
      }
    }

    //Rule 2: relabel national cities based on proximity
    for (size_t i = 0; i < vertex_table.size(); i++) {
      if (vertex_table[i].getType() == NATIONAL) {
        //find all national cities withing 150 miles
        vector<int> national_indicies;
        for (size_t j = 0; j < vertex_table.size(); j++) {
          if (i != j && vertex_table[j].getType() == NATIONAL && dist_table.getValue(i,j) < 150.0) {
            national_indicies.push_back(j);
          }
        }

        //find national city with largest pop
        int largest_pop_index = i;
        for (int index : national_indicies) {
          if (vertex_table[index].getPopulation() > vertex_table[largest_pop_index].getPopulation()) {
            largest_pop_index = index;
          }
        }

        //relabel other national cities as regional
        for (int index : national_indicies) {
          if (index != largest_pop_index) {
            vertex_table[index].setType(REGIONAL);
          }
        }
      }
    }

    //Rule 3: relabel regional cities based on proximity
    for (size_t i = 0; i < vertex_table.size(); i++) {
      if (vertex_table[i].getType() == REGIONAL) {
        //find regional cities within 100 miles
        vector<int> regional_indicies;
        for (size_t j = 0; j < vertex_table.size(); j++) {
          if (i != j && vertex_table[j].getType() == REGIONAL && dist_table.getValue(i,j) < 100.0) {
            regional_indicies.push_back(j);
          }
        }

        //find largest regional city
        int largest_pop_index = i;
        for (int index : regional_indicies) {
          if (vertex_table[index].getPopulation() > vertex_table[largest_pop_index].getPopulation()) {
            largest_pop_index = index;
          }
        }

        //relabel other regional cities as local
        for (int index : regional_indicies) {
          if (index != largest_pop_index) {
            vertex_table[index].setType(LOCAL);
          }
        }
      }
    }
}

//calculate distance from given equation
float calcDistance(const city& city1, const city& city2) {
    float lat1 = city1.getLatitude() * DEG2RAD;
    float lon1 = city1.getLongitude() * DEG2RAD;
    float lat2 = city2.getLatitude() * DEG2RAD;
    float lon2 = city2.getLongitude() * DEG2RAD;

    float dlon = lon2 - lon1;
    float dlat = lat2 - lat1;

    float a = pow(sin(dlat/2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlon/2), 2);
    float c = 2 * asin(sqrt(a));

    return earthradius * c;
}

void create_dist_table(vector<city>& vertex_table, matrix &dist_table) {
    int num_cities = vertex_table.size();
    //iterate over pairs of cities
    for (int i = 0; i < num_cities; i++) {
      for (int j = 0; j < num_cities; j++) {
        float distance = calcDistance(vertex_table[i], vertex_table[j]);
        //round and store distance in dist_table
        distance = 5.0 * round(distance / 5.0);
        dist_table.setValue(i,j, distance);
      }
    }
}
void create_time_table(vector<city> &vertex_table, matrix &dist_table, matrix &time_table) {

    int num_cities = vertex_table.size();

    //conversion rates
    const float plane_speed = 520.0;
    const float truck_speed = 65.0;

    //iterate over pairs of cities
    for (int i = 0; i < num_cities; i++) {
      for (int j = 0; j < num_cities; j++) {  
        float distance = dist_table.getValue(i,j);

        //calc travel time
        float travel_time;
        if (vertex_table[i].getType() == NATIONAL && vertex_table[j].getType() == NATIONAL) {
          //national cities by plane
          travel_time = distance / plane_speed;
        } else {
          //other types by truck
          travel_time = distance / truck_speed;
        }

        //add times to time table
        time_table.setValue(i, j, travel_time);
      }
    }
}

bool compare_dist(const pair<float, int> &a, const pair<float, int> &b) {
  return a.first < b.first;
}

void create_edge_table(vector<city> &vertex_table, matrix &edge_table, matrix &dist_table) { 

    int num_cities = vertex_table.size();

    //Rule 1: National to each other but not themselves
    for (int i = 0; i < num_cities; i++) {
      if (vertex_table[i].getType() == NATIONAL) {
        for (int j = 0; j < num_cities; j++) {
          if (i != j && vertex_table[j].getType() == NATIONAL) {
            float distance = dist_table.getValue(i,j);
            
            //connect cities bidirectionall and store distance in edge table
            edge_table.setValue(i, j, distance);
            edge_table.setValue(j, i, distance);
          }
        }
      }
    }

    //Rule 2: Connect regional cities to 3 nearest non local cities
    for (int i = 0; i < num_cities; i++) {
      if (vertex_table[i].getType() == REGIONAL) {
        vector<pair<float, int>> distance_map;//store distance and indicies

        //create map to store distances from current city to other cities
        for (int j = 0; j < num_cities; j++) {
          if (i!=j && vertex_table[j].getType() != LOCAL) { //exclude local cities
            float distance = dist_table.getValue(i, j);
            distance_map.push_back(make_pair(distance, j));
          }
        }

        //sort first 3 by lowest distance
        partial_sort(distance_map.begin(), distance_map.begin() + 3, distance_map.end(), compare_dist);
        
        int count = 0;
        //connect those 3 cities
        for (auto iT = distance_map.begin(); count < 3 && iT != distance_map.end(); iT++) {
          int index = iT->second;
          float distance = iT->first;

          //connect cities both directions and store distance in edge table
          edge_table.setValue(i, index, distance);
          edge_table.setValue(index, i, distance);
          count++;
        }
      }
    }

    //Rule 3: Connect LOCAL cities to 5 nearest non local and any local under 50 miles
    for (int i = 0; i < num_cities; i++) {
      if (vertex_table[i].getType() == LOCAL) {
        vector<pair<float, int>> distance_map; 
        for (int j = 0; j < num_cities; j++) {
          if (i != j && vertex_table[j].getType() != LOCAL) {
            float distance = dist_table.getValue(i , j);
            distance_map.push_back(make_pair(distance, j));
          }
        }

        //sort 5 closest
        partial_sort(distance_map.begin(), distance_map.begin() + 5, distance_map.end(), compare_dist);

          //connect to 5 non local cities
          int count = 0;
          
          for (auto iT = distance_map.begin(); count < 5 && iT != distance_map.end(); iT++) {
            int index = iT->second;
            float distance = iT->first;

            edge_table.setValue(i, index, distance);
            edge_table.setValue(index, i, distance);
            count++;
          }

          //connect any local cities less than 50 miles

          for (int j = 0; j < num_cities; j++) {
            if (i != j && vertex_table[j].getType() == LOCAL && dist_table.getValue(i,j) < 50.0) {
              float distance = dist_table.getValue(i,j);
              //connect both directions
              edge_table.setValue(i, j, distance);
              edge_table.setValue(j, i, distance);

            }
          }
        
      }
    }
}

void write_vertex_table(vector<city>& vertex_table) { 

    ofstream file("vertex_table.txt");
    if (!file.is_open()) {
      cerr << "unable to open vertex table for writing" << endl;
      return;
    }

    for (int i = 0; i < (int)vertex_table.size(); i++) {
      file << setw(4) << right << i << "  " << vertex_table[i] << endl;
    }

    file.close();
}

void write_dist_table(vector<city> &vertex_table,matrix &dist_table) { 
    
    ofstream file("dist_table.txt");
    if (!file.is_open()) {
      cerr << "unable to open dist table for writing" << endl;
      return;
    }

    for (int i = 0; i < dist_table.getSize(); i++) {
      for (int j = 0; j < i; j++) {
        if (i != j) {
                float distance = dist_table.getValue(i, j);
                file << setw(4) << i  << "  " << setfill('.') << setw(19) << left << vertex_table[i].getName() << " to "
                        << setw(19) << left << vertex_table[j].getName() << setfill(' ') << setw(8) << right
                        << fixed << setprecision(1) << distance << " miles" << endl;
        }
      }
      if (i != 0) {
        file << "\n";
      }
    }
  

    file.close();
}

void write_time_table(vector<city> &vertex_table,matrix &time_table) { 
    ofstream file("time_table.txt");
    if (!file.is_open()) {
      cerr << "unable to open time table for writing" << endl;
      return;
    }

    for (int i = 0; i < time_table.getSize(); ++i) {
        for (int j = 0; j < i; ++j) {
            float time = time_table.getValue(i, j);
            file << setw(4) << i << "  " << setfill('.') <<setw(19) << left << vertex_table[i].getName() << " to "
                 << setw(19) << left << vertex_table[j].getName()  << setfill(' ') << setw(8) << right
                 << fixed << setprecision(1) << time << " hours" << endl;
        }
        if (i != 0) {
          file << "\n";
        }
    }
    file.close();
}

void write_edge_table(vector<city> &vertex_table, matrix &edge_table,matrix &dist_table, matrix &time_table) { 
    ofstream file("edge_table.txt");
    if (!file.is_open()) {
      cerr << "unable to open edge table for writing" << endl;
      return;
    }
    
    int num_cities = edge_table.getSize();
  
    for (int i = 0; i < num_cities; i++) {
      file << setw(4) << i << " " << vertex_table[i].getName() << endl;
      for (int j = 0; j < num_cities; j++) {
        if (edge_table.getValue(i, j) > 0) {
             file << "  " << right << setw(4) << j << " " << setfill('.') << vertex_table[j].getName();
             file << setw(city::w - vertex_table[j].getName().length() + 5);
             file << "  " << setfill(' ')
                  << fixed << setprecision(1) << setw(6) << right << dist_table.getValue(i, j) << " miles"
                  << fixed << setprecision(1) << setw(5) << right << time_table.getValue(i, j) << " hours" << endl;
        }
      }
      if (i != num_cities - 1) {
        file << "\n";
      }
    }
    

    file.close();

}

//return city index for city name for dijkstra
int city_index(vector<city> &vertex_table, const string &city_name) {
  for (size_t i = 0; i < vertex_table.size(); i++) {
    if (vertex_table[i].getName() == city_name) {
      return i;
    }
    
  }
  //if not found, return -1
  return -1;
}


//implement code from graph handout 3 to show dijkjstra route for time and distance
void dijkstra_route(int city_from, int city_to, vector<city> &vertex_table, matrix &edge_table,
                    string mode, matrix &dist_table, matrix &time_table) { 
    
    int num_cities = vertex_table.size();

    vector<color_t> v_color(num_cities, WHITE);  //vertex color
    vector<float> distance(num_cities, numeric_limits<float>::max());  //shortest distance
    vector<int> parent(num_cities, -1);  //parent of each vertex in the shortest path

    // Initialize the starting vertex
    distance[city_from] = 0;

    
    for (int i = 0; i < num_cities; ++i) {
        //find vertex with min distance
        int cur_v = -1;
        float min_dist = numeric_limits<float>::max();
        for (int j = 0; j < num_cities; ++j) {
            if (v_color[j] == WHITE && distance[j] < min_dist) {
                cur_v = j;
                min_dist = distance[j];
            }
        }

        //if no unvisited vertex found break
        if (cur_v == -1) break;

        //mark current vertex as black
        v_color[cur_v] = BLACK;

        //update distance to adj vertices
        for (int adj_v = 0; adj_v < num_cities; ++adj_v) {
            float weight;
            if (mode == "-time") {
                weight = time_table.getValue(cur_v, adj_v);  
            } else if (mode == "-dist") {
                weight = dist_table.getValue(cur_v, adj_v);  
            }

          if (edge_table.getValue(cur_v,adj_v) != 0) {  
            if (v_color[adj_v] == WHITE && distance[cur_v] + weight < distance[adj_v]) {
                distance[adj_v] = distance[cur_v] + weight;
                parent[adj_v] = cur_v;
            }
          }
        }
    }

    //output the shortest path from city_from to city_to
    int current = city_to;
    vector<int> d_path;
    float total_distance = 0.0;
    float total_time = 0.0;

    
    while (current != -1) {
        d_path.push_back(current);
        current = parent[current];
    }

    
    
    // Print the path in reverse order
    bool firstit = true;
    for (int i = d_path.size() - 1; i >= 0; --i) {
        
        int city_index = d_path[i];

        //incremental distance and time
        float inc_dist = 0.0;
        float inc_time = 0.0;

        
        if (!firstit) { //if not first iteration, calc incremental distance and time for all iterations except first one
        int prev_city_index = d_path[i + 1];
        inc_dist = dist_table.getValue(prev_city_index, city_index);
        inc_time = time_table.getValue(prev_city_index, city_index);
        }

        //calc total times and distance
        total_distance += inc_dist;
        total_time += inc_time;

        //output path with incremental dist and times
        cout << setw(8) << fixed << setfill(' ') << setprecision(2) << right << total_distance << " miles " 
             << setw(5) << fixed << setprecision(2) << right << total_time <<  " hours "
             << setw(4) << city_index << " " << vertex_table[city_index].getName();

             int dots = city::w - vertex_table[city_index].getName().length() + 3;
        for (int j = 0; j < dots; j++) {
            cout << ".";
        }
        
        cout << "  " << setw(8) << right << vertex_table[city_index].getPopulation();

          if (!firstit) { //print incremental distance and time for all iterations except first one
            cout << setw(8) << fixed << setprecision(2) << right << inc_dist << " miles "
            << setw(5) << fixed << setprecision(2) << right << inc_time << " hours";
        }

        
          cout << endl;
     
          if (firstit) {
            firstit = false;
          }
    }
    cout << "\n";
} 

int main(int argc, char *argv[]) {
  //parse cmd line options
  if (argc != 3) {
    cerr << "usage: ./Prog6 -info|dist|time [-seed=N] cities.txt" << endl;
    return 1;
  }

  string arg = argv[1];
  string filename = argv[argc-1];


  if (arg!="-info" && arg!="-dist" && arg!="-time") {
    cerr << "usage: ./Prog6 -info|dist|time [-seed=N] cities.txt" << endl;
    return 1;
  }

  ifstream infile(filename);

  //read data and create vertex table
  vector<city> vertex_table;
  create_vertex_table(infile, vertex_table);
  infile.close();

  //create and fill dist table
  matrix dist_table(vertex_table.size());
  create_dist_table(vertex_table, dist_table);

  //update vertex table based on pop and proximity rules
  update_vertex_table(vertex_table, dist_table);

  //create and fill time table
  matrix time_table(vertex_table.size());
  create_time_table(vertex_table, dist_table, time_table);

  //create edge table
  matrix edge_table(vertex_table.size());
  create_edge_table(vertex_table, edge_table, dist_table);

  //for info write tables
  if (arg == "-info") {
    write_vertex_table(vertex_table);
    write_dist_table(vertex_table, dist_table);
    write_time_table(vertex_table, time_table);
    write_edge_table(vertex_table, edge_table,dist_table, time_table);
    return 0;
  }

  //for dijkstra, ask for city from and city to and output shortest path
  bool dijkstra = true;
  while (dijkstra) {
    string city_from, city_to;
    cout << "Enter> ";
    cin >> city_from >> city_to;

    //eof(ctrl d) break
    if (cin.eof()) {
      cout << endl;
      break;
    }

    int city_from_i = city_index(vertex_table, city_from);
    int city_to_i = city_index(vertex_table, city_to);

    if (city_from_i == -1 || city_to_i == -1) {
            cerr << "Invalid_city: prefix not found" << endl;
            continue; // Retry input
    }

    string mode; 

    //set mode for dijkstra
    if (arg == "-dist") {
      mode = "-dist";
    } else
    if (arg == "-time") {
      mode = "-time";
    }


    dijkstra_route(city_from_i, city_to_i, vertex_table, edge_table, mode, dist_table, time_table);
  }
  return 0;
}
