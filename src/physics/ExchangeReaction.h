/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __EXCHANGE_REACTION_H__
#define __EXCHANGE_REACTION_H__

#include <random>
#include <tuple>
#include <unordered_map>
#include <sstream>
#include <string>

#include <easyloggingcpp/easylogging++.h>

#include <physics/Element.h>
#include <physics/SimulationBox.h>

class ExchangeReaction {
private:

    //! Logger name.
    const char *logName_;

    //! random number generator
    //std::default_random_engine random_;
    std::mt19937 random_;

    //! Exchange probabilities
    std::unordered_map<std::string, double> reactions_;

    //! Generates unique ID for specified element IDs.
    inline std::string getId_(const std::string material1,
                              const std::string material2) const;

    //! Extracts element ID from unique ID for reaction
    inline std::string getMaterial_(const std::string id, bool first);

public:
    //! Constructor
    ExchangeReaction (const char *logName="ExchangeReaction");

    //! Destructor
    ~ExchangeReaction ();

    //! Adds an entry for an exchange reaction with a given probability
    void addReaction (const std::string material1, const std::string material2,
                      const double probability);
    //! Adds an entry for an exchange reaction with a given probability
    void addReaction (const std::string reactionDefinition);

    //! Replaces elements in the specified layer according to the defined reactions
    void replaceAnions (SimulationBox *simbox, const indexType layerID,
                        const indexType dimension,
                        const std::string original, const std::string replacement);

    //! Returns information about stored reactions
    std::string str ();
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
