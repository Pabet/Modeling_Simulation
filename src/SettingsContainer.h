//
// Created by pabet on 12.12.18.
//

#ifndef PSEMOLDYN_GROUPG_SETTINGSGENERATOR_H
#define PSEMOLDYN_GROUPG_SETTINGSGENERATOR_H

#include "utils/Vector.h"
#include <list>

class SettingsContainer {
public:
    double ent_time;
    double delta_t;
    double factor;
    utils::Vector<double, 3> domain_size;
    double rcutoff;
    int x1_boundary_condition;
    int x2_boundary_condition;
    int y1_boundary_condition;
    int y2_boundary_condition;
    int z1_boundary_condition;
    int z2_boundary_condition;
    bool brownian_motion;
    double initial_temperature;
    double nthermostat;
    double target_temperature;
    double temperature_difference;
    double gravitation;
    bool writeCheckpoint;
    bool readCheckpoint;
    double writeCheckpointTime;
    int parallelisation_method;
    int force_calculation_method;
    double r0;
    double k;
    double fz_up;
    double rl;


public:
    SettingsContainer(double &et, double &dt, double &f, utils::Vector<double, 3> &ds,
                      double &rcutoff, bool &brownian_motion, double &initial_temperature,
                      double &nthermostat, double &target_temperature, double &temperature_difference, double &gravitation, bool &write_checkpoint, bool &read_checkpoint, double &write_checkpoint_time,
                      int &x1, int &x2, int &y1, int &y2, int &z1, int &z2, double &r0v, double &kv, double &fz_upv, int &pm, int &fcm, double &rl);

    SettingsContainer() = default;

    double getEnt_time() const;

    double getFactor() const;

    const utils::Vector<double, 3> &getDomain_size() const;

    double getRcutoff() const;

    bool isBoundary_condition() const;

    bool isBrownian_motion() const;

    double getInitial_temperature() const;

    double getNthermostat() const;

    double getTarget_temperature() const;

    double getTemperature_difference() const;

    double getGravitation() const;

    bool getWriteCheckpoint () const;

    bool getReadCheckpoint () const;

    double getWriteCheckpointTime () const;
};


#endif //PSEMOLDYN_GROUPG_SETTINGSGENERATOR_H
