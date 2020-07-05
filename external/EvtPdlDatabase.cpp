#include "EvtPdlDatabase.h"
#include "StringTools.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <tuple>
#include <stdexcept>
#include <string_view>

using namespace std;

namespace sct::comm {

string EvtPdlDatabase::evtPdlFile("../resource/evt.pdl");
hash<string> hasher;
optional<const EvtPdlDatabase> EvtPdlDatabase::instance;

tuple<EvtPdlItem, string, int32_t> lineToPdlItem(string_view sw) {
    auto parts = ST::SplitAndStrip(sw, " ");
    if (parts.size() != 12) {
        cerr << "EvtPdlDatabase: wrong file format.\n"
             << " line: " << sw << endl;
        throw std::runtime_error("EvtPdlDatabase:  wrong file format. " + string(sw));
    }
    string name{parts[3]};
    size_t hname = hasher(name);
    try {
        int32_t id = stoi(string(parts[4]));
        double mass = stod(string(parts[5]));
        double width = stod(string(parts[6]));
        int8_t charge = stoi(string(parts[8]));
        uint8_t spin = stoi(string(parts[9]));
        double lifetime = stod(string(parts[10]));

        EvtPdlItem item{mass, width, lifetime, hname, charge, spin};
        return tie(item, name, id);
    } catch (invalid_argument& ex) {
        cerr << "EvtPdlDatabase: wrong file format.\n"
             << " invalid argument: " << ex.what() << endl;
        throw ex;
    } catch (out_of_range& ex) {
        cerr << "EvtPdlDatabase: wrong file format.\n"
             << " out_of_range: " << ex.what() << endl;
        throw ex;
    }
}

EvtPdlDatabase::EvtPdlDatabase(const std::string& fname) {
    ifstream ifile(fname, ifstream::in);
    if (!ifile.good()) {
        cerr << "EvtPdlDatabase: Can't open file " << fname << endl;
        throw std::runtime_error("EvtPdlDatabase: Can't open file " + fname);
    }

    string line;
    while(!ifile.eof()) {
        getline(ifile, line);
        if (line.substr(0, 3) == "add") {
            auto [item, name, id] = lineToPdlItem(line);
            hname_to_id[item.hname] = id;
            hashed_names[item.hname] = move(name);
            id_to_item[id] = move(item);
        }
    }
}

const EvtPdlDatabase& EvtPdlDatabase::Instance() {
    if (!instance.has_value()) {
        instance.emplace(EvtPdlDatabase(evtPdlFile));
    }
    return instance.value();
}

bool EvtPdlDatabase::isKnownId(int32_t id) const {
    return id_to_item.find(id) != id_to_item.end();
}

bool EvtPdlDatabase::isKnownName(const std::string& name) const {
    return hname_to_id.find(hasher(name)) != hname_to_id.end();
}

const string& EvtPdlDatabase::name(int32_t id) const {
    return hashed_names.at(id_to_item.at(id).hname);
}

int32_t EvtPdlDatabase::id(const string& name) const {
    return hname_to_id.at(hasher(name));
}

int32_t EvtPdlDatabase::antiParticleId(int32_t id) const {
    if (auto item = id_to_item.find(-id); item != id_to_item.end()) {
        return -id;
    }
    return id;
}

string EvtPdlDatabase::antiParticleName(const string& name) const {
    return hashed_names.at(
        id_to_item.at(antiParticleId(id(name))).hname
    );
}

double EvtPdlDatabase::mass(int32_t id) const {
    return id_to_item.at(id).mass;
}

double EvtPdlDatabase::mass(const string& name) const {
    return mass(id(name));
}

double EvtPdlDatabase::width(int32_t id) const {
    return id_to_item.at(id).width;
}

double EvtPdlDatabase::width(const string& name) const {
    return width(id(name));
}

double EvtPdlDatabase::spin(int32_t id) const {
    return 0.5 * id_to_item.at(id).spin;
}

double EvtPdlDatabase::spin(const string& name) const {
    return spin(id(name));
}

double EvtPdlDatabase::lifetime(int32_t id) const {
    return id_to_item.at(id).lifetime;
}

double EvtPdlDatabase::lifetime(const string& name) const {
    return lifetime(id(name));
}

double EvtPdlDatabase::charge(int32_t id) const {
    return id_to_item.at(id).charge / 3.;
}

double EvtPdlDatabase::charge(const string& name) const {
    return charge(id(name));
}

}  // namespace sct::comm
