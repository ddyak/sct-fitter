/**************************************************************************
 * AURORA (SCT Framework)                                                *
 * Copyright(C) 2019 - SCT Collaboration                                  *
 *                                                                        *
 * Author: The SCT Collaboration                                          *
 * Contributors: Vitaly Vorobyev                                          *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/
#pragma once

#include <unordered_map>
#include <string>
#include <cstdint>
#include <functional>
#include <optional>

namespace sct::comm {

struct EvtPdlItem {
    double mass;
    double width;
    double lifetime;
    size_t hname;
    int8_t charge;  // 3*charge
    uint8_t spin;   // 2*spin
};

class EvtPdlDatabase {
    std::unordered_map<size_t, std::string> hashed_names;
    std::unordered_map<size_t, int32_t> hname_to_id;
    std::unordered_map<int32_t, EvtPdlItem> id_to_item;

    static std::optional<const EvtPdlDatabase> instance;

    EvtPdlDatabase(const std::string& fname);  // singleton

 public:
    static std::string evtPdlFile;
    static const EvtPdlDatabase& Instance();

    bool isKnownId(int32_t id) const;
    bool isKnownName(const std::string& name) const;

    const std::string& name(int32_t id) const;
    int32_t id(const std::string& name) const;
    
    int32_t antiParticleId(int32_t id) const;
    std::string antiParticleName(const std::string& name) const;

    double mass(int32_t id) const;
    double mass(const std::string& name) const;

    double width(int32_t id) const;
    double width(const std::string& name) const;

    double spin(int32_t id) const;
    double spin(const std::string& name) const;

    double lifetime(int32_t id) const;
    double lifetime(const std::string& name) const;

    double charge(int32_t id) const;
    double charge(const std::string& name) const;
};

}  // namespace sct::comm
