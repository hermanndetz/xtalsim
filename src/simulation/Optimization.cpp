/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "Optimization.h"

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
Optimization::Optimization(const char *logName): logName_(logName)
{
    random_.seed(std::random_device()());
}

//------------------------------------------------------------------------------

Optimization::~Optimization() {}

//------------------------------------------------------------------------------

//! Adds a function to the list of callable optimizations. The given probability
//! is used to determine how often compared to other optimizations, this one is
//! called.
//! \param name Name of registered action. Is later used for example in the
//! statistic.
//! \param fp Function pointer to the function which shall be called.
//! \param probability Probability with which the function shall be called.
void Optimization::registerAction(const std::string name,
				  const OptimizationFunction fp,
				  const double probability)
{
    CLOG(TRACE, logName_) << "beginning to register Action";
    
    actions_.push_back(fp);
    statistics_.push_back(OptimizationStatistic(name));
    probabilities_.push_back(probability);
    probabilitiesSum_.push_back(probability);
    checkProbabilities();

    CLOG(TRACE, logName_) << "Action sucessfully registered";
}

//------------------------------------------------------------------------------

//! Adapts probabilities such that the sum of all is equal to 1. In detail the
//! sum of all probablities is computed. In the next step each single one of
//! them starting at the first is multiplied and added to a sum. This sum is in
//! each step written to a separate vector. It therefore contains the sum of all
//! probabilieties up to this point, whereat it is guaranteed that the last
//! entry is 1.
void Optimization::checkProbabilities(void)
{
    double sum=0, factor;

    CLOG(TRACE, logName_) << "beginning to check probabilities";
    
    for(auto el: probabilities_)
	sum+= el;
    
    factor = 1/sum;

    CLOG(DEBUG, logName_) << "got probability sum of " << sum;
    CLOG(DEBUG, logName_) << "applying factor " << factor <<
	" to single probabilities";
    
    sum = 0;

    std::transform(probabilities_.begin(), probabilities_.end(),
		   probabilitiesSum_.begin(),
		   [&sum,factor] (double prob) {
		       sum += (prob*factor);
		       return sum;
		   } );
}

//------------------------------------------------------------------------------

//! Pick a random action from the available ones and return the appropriate
//! index. More specific a random value between 0 and 1 is picked. Then the
//! vector probabilitiesSum_, which contains the scaled sum of the
//! probabilities, is traversed from the start. If a value > the random choice
//! is detected the random choice was achieved.
//! \return Randomly chosen vector index.
inline uint16_t Optimization::getRandomAction(void)
{

    
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double choice = distribution(random_);

    CLOG(TRACE, logName_) << "get random action " << choice;
    
    // std::for_each(probabilities_.begin(), probabilities_.end(),
    // 	    [this] (double prob) {logger_->info("got probability {}", prob);});

    // std::for_each(probabilitiesSum_.begin(), probabilitiesSum_.end(),
    // 	    [this] (double prob) {logger_->info("got prob sum {}", prob);});

    for (uint16_t i=0; i<probabilitiesSum_.size(); i++){
	if (probabilitiesSum_[i] > choice)
	    return(i);
    }

    return(actions_.size()-1);
}

//------------------------------------------------------------------------------

//! Runs the optimizations until a certain precision is reached. For this task
//! the energy gains are checked and, if necessary, the parameters are adapted
//! approriately to reach the desired goal.
//! \param simbox Reference to the simulation box that is optimized.
//! \param range Range of the simulation box that is optimized.
//! \param tersoff Parameters necessary to compute the energy.
//! \param parameter Optimization parameters such as run count or maximal
//! displacement.
//! \param journalPreamble Path preamble that shall be used to store generated
//! journals.
void Optimization::runDynamic(SimulationBox &simbox,
		       const Range3D<indexType> &range,
		       const TersoffPotential &tersoff,
		       OptimizationParameter parameter,
		       const std::string journalPreamble)
{

    int action, counter=0;
    XmlHandler xmlHandler;
    double startEnergy, refEnergy;
    std::string journalFileName = journalPreamble + "optimizationEnergy.xml";
    Journal<double> energyJournal("optStepEnergy","Energy after each"
				  "optimization step");
    std::ofstream energyFile;

    CLOG(TRACE, logName_) << "Starting optimized optimization";
    CLOG(INFO, logName_) << "Checking each " << parameter.checkCount <<
	" steps if energy dropped by factor " << parameter.energyDropFactor <<
	". Parameters are reduced with factor " << parameter.reductionFactor <<
	" until minimum of MMC(" << parameter.mmcMinDisplacement <<
	") or Scaling (" << parameter.minScaling <<") was reached or energy " <<
	"gain in first runs was below " << parameter.minEnergy;

    // is updated with current energy in real time
    energyFile.open(journalPreamble + "energies",
		    std::ios::out | std::ios::app);
    energyFile.precision(20);
    energyFile << "####" << std::endl;
    
    try{
	xmlHandler.load(journalFileName);
	xmlHandler.get(energyJournal,"optStepEnergy");
    }catch(XmlException e){
	//file not available, do nothing
    }
    
    double energy = simbox.getEnergy(tersoff);
    CLOG(DEBUG, logName_) << "energy after start: " << energy;
    energyJournal.add(energy);

    counter =0;
    refEnergy=0;
    while (true){

	if (counter ==0)
	    startEnergy = simbox.getEnergy(tersoff);
	
	action=getRandomAction();
	counter++;
	
	(simbox.*( actions_[action] ))(tersoff, parameter,
				       statistics_[action], range);

	double energy = simbox.getEnergy(tersoff);
	CLOG(DEBUG, logName_) << "energy after step: " << counter << ": "
			      << energy;
	energyJournal.add(energy);
	energyFile << energy << std::endl;
	
	if (counter != parameter.checkCount) continue;

	CLOG(TRACE, logName_) << "checking energy gain";
	double energyGain = startEnergy - simbox.getEnergy(tersoff);
	counter=0;
	    
	if (refEnergy > 0){
	    if (energyGain < refEnergy){
		CLOG(DEBUG, logName_) << "energy gain was below limit";
		CLOG(INFO, logName_) << "reducing parameters after " <<
		    counter << " steps";
		CLOG(INFO, logName_) << "energy gain was " << energyGain <<
		    "compared to reference " << refEnergy;

		parameter.mmcMaxDisplacement /= parameter.reductionFactor;
		parameter.maxScaling /= parameter.reductionFactor;
		CLOG(DEBUG, logName_) << "using max displacement of " <<
		    parameter.mmcMaxDisplacement << " and max Scaling of " <<
		    parameter.maxScaling;
		this->printStatistic();
		this->clearStatistic();

		if (parameter.mmcMaxDisplacement < parameter.mmcMinDisplacement)
		    break;
		    
		if (parameter.maxScaling < parameter.minScaling)
		    break;
		    
		refEnergy = 0;

	    }
	    else{
		CLOG(DEBUG, logName_) << "energy gain was above limit";
		CLOG(INFO, logName_) << "energy gain was " << energyGain <<
		    " compared to reference " << refEnergy;
	    }
	}
	else {
	    refEnergy = energyGain/parameter.energyDropFactor;
	    CLOG(INFO, logName_) << "reference energy gain set to "<< refEnergy;

	    if (refEnergy < parameter.minEnergy){
		CLOG(INFO, logName_) << "energy gain achieved is below " <<
		    "minimal energy";
		break;
	    }
	}

    }

    energyFile.close();
    xmlHandler.set(energyJournal);
    xmlHandler.save(journalFileName);
    CLOG(TRACE, logName_) << "successfully finished optimized optimization";
} 

//------------------------------------------------------------------------------
    
//! Runs the optimizations.
//! \param simbox Reference to the simulation box that is optimized.
//! \param range Range of the simulation box that is optimized.
//! \param tersoff Parameters necessary to compute the energy.
//! \param parameter Optimization parameters such as run count or maximal
//! displacement.
//! \param journalPreamble Path preamble that shall be used to store generated
//! journals.
void Optimization::runStatic(SimulationBox &simbox,
                             const Range3D<indexType> &range,
                             const TersoffPotential &tersoff,
                             const OptimizationParameter &parameter,
                             const std::string journalPreamble)
{

    int reportCount, reportCountStep, action;
    Journal<double> energyJournal("optStepEnergy","Energy after each"
                                  "optimization step");

    Journal<UDTuple> eJournal("eJournal","Energy after each optimization step");

    XmlHandler xmlHandler;
    XmlHandler xmlHandler2;

    std::string journalFileName = journalPreamble + "optimizationEnergy.xml";
    std::string journalFileName2 = journalPreamble + "optEnergy.xml";
    reportCount = reportCountStep = parameter.runCount/10;
    bool doReport = (reportCountStep > 0);

    CLOG(TRACE, logName_) << "Starting optimization";
    CLOG(DEBUG, logName_) << "Overall " << parameter.runCount <<" runs are done";

    if (parameter.runCount < 1){
        CLOG(WARNING, logName_) << "Run count below one! Doing nothing!";
        return;
    }

    if (journalPreamble != ""){
        try{
            xmlHandler.load(journalFileName);
            xmlHandler.get(energyJournal,"optStepEnergy");
            xmlHandler2.load(journalFileName2);
            xmlHandler2.get(eJournal,"eJournal");
        }catch(XmlException e){
            //file not available, do nothing
        }
    }

    double energy = simbox.getEnergy(tersoff);
    CLOG(DEBUG, logName_) << "Energy after start: " << energy;
    energyJournal.add(energy);
    eJournal.add(std::make_tuple(0, energy/simbox.getAtomCount()));
    
    for (int i=0; i<parameter.runCount; i++){

        action=getRandomAction();
	
        (simbox.*( actions_[action] ))(tersoff, parameter,
                                       statistics_[action], range);
        //actions_[getRandomAction()](range, tersoff, value); // std::function

        if (doReport &&(reportCount == i)){
            CLOG(INFO, logName_) << "finished optimization # " << i;
            reportCount += reportCountStep;
        }

        double energy = simbox.getEnergy(tersoff);
        CLOG(DEBUG, logName_) << "energy after step " << i << ": " << energy;
        energyJournal.add(energy);
        eJournal.add(std::make_tuple(i, energy/simbox.getAtomCount()));
    }

    if (journalPreamble != ""){
        xmlHandler.set(energyJournal);
        xmlHandler.save(journalFileName);

        xmlHandler2.set(eJournal);
        xmlHandler2.save(journalFileName2);
    }
} 

//------------------------------------------------------------------------------

//! Prints the name, probability, number of accepted, rejected and overall
//! tries. 
void Optimization::printStatistic(void) const
{

    CLOG(INFO, logName_) << "Statistic of optimization actions:";
    
    for (unsigned int i=0; i<statistics_.size(); i++){
	auto el = statistics_[i];
	CLOG(INFO, logName_) << "Action name: " << el.name;
	
	if (i==0){
	    CLOG(INFO, logName_) << "Probability: " << probabilitiesSum_[0];
	}
	else{
	    CLOG(INFO, logName_) << "Probability: " << probabilitiesSum_[i]-
		probabilitiesSum_[i-1];
	}
	
	CLOG(INFO, logName_) << "Accept count: " << el.acceptCount;
	CLOG(INFO, logName_) << "Reject count: " << el.rejectCount;
	CLOG(INFO, logName_) << "Overall tries: " <<
	    el.rejectCount+el.acceptCount;	
    }
}

//------------------------------------------------------------------------------

//! Set the count of accepted and rejected moves back to 0.
void Optimization::clearStatistic(void)
{

    CLOG(TRACE, logName_) << "Clearing optimization statistic.";
    
    for (unsigned int i; i<statistics_.size(); i++){
	statistics_[i].acceptCount=0;
	statistics_[i].rejectCount=0;
    }

    CLOG(TRACE, logName_) << "Statistic cleared.";
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
