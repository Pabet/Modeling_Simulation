//
// Created by pabet on 12.12.18.
//

#include "SettingsContainer.h"
#include "utils/Vector.h"
#include <list>


SettingsContainer::SettingsContainer (double &et, double &dt, double &f, utils::Vector<double, 3> &ds, double &rc,
                                      bool &bm, double &it, double &nt, double &tt, double &td, double &g, bool &wc,
                                      bool &rch, double &wct, int &x1, int &x2, int &y1, int &y2, int &z1, int &z2,
                                      double &r0v, double &kv, double &fz_upv, int &pm, int &fcm, double &rlv)
                                     {
    ent_time = et;
    delta_t = dt;
    factor = f;
    domain_size.operator = (ds);
    rcutoff = rc;
    x1_boundary_condition = x1;
    x2_boundary_condition = x2;
    y1_boundary_condition = y1;
    y2_boundary_condition = y2;
    z1_boundary_condition = z1;
    z2_boundary_condition = z2;
    brownian_motion = bm;
    initial_temperature = it;
    nthermostat = nt;
    target_temperature = tt;
    temperature_difference = td;
    gravitation = g;
    writeCheckpoint = wc;
    readCheckpoint = rch;
    writeCheckpointTime = wct;
    r0 = r0v;
    k = kv;
    fz_up = fz_upv;
    parallelisation_method = pm;
    force_calculation_method = fcm;
    rl = rlv;
}

double SettingsContainer::getEnt_time() const {
    return ent_time;
}

double SettingsContainer::getFactor() const {
    return factor;
}

const utils::Vector<double, 3> &SettingsContainer::getDomain_size() const {
    return domain_size;
}

double SettingsContainer::getRcutoff() const {
    return rcutoff;
}

bool SettingsContainer::isBrownian_motion() const {
    return brownian_motion;
}

double SettingsContainer::getInitial_temperature() const {
    return initial_temperature;
}

double SettingsContainer::getNthermostat() const {
    return nthermostat;
}

double SettingsContainer::getTarget_temperature() const {
    return target_temperature;
}

double SettingsContainer::getTemperature_difference() const {
    return temperature_difference;
}

double SettingsContainer::getGravitation() const {
    return gravitation;
}

bool SettingsContainer::getWriteCheckpoint () const {
    return writeCheckpoint;
}

bool SettingsContainer::getReadCheckpoint() const {
    return readCheckpoint;
}

double SettingsContainer::getWriteCheckpointTime() const {
    return writeCheckpointTime;
}