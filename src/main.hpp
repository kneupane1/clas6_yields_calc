/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "color.hpp"
#include "histogram.hpp"
#include "reaction.hpp"

size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists, int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();
  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  auto data = std::make_shared<Branches6>(_chain);

  // Total number of events "Processed"
  size_t total = 0;
  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);
    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    // Quick electron cuts
    bool elec = true;
    elec &= (data->q(0) == NEGATIVE);
    elec &= (data->ec(0) > 0);
    elec &= (data->gpart() > 0);
    elec &= (data->gpart() < 7);
    elec &= (data->cc(0) > 0);
    elec &= (data->stat(0) > 0);
    elec &= (data->sc(0) > 0);
    elec &= (data->dc(0) > 0);
    elec &= (data->dc_stat(0) > 0);
    if (!elec) continue;

    // If we pass electron cuts the event is processed
    total++;

    // Make a reaction class from the data given
    auto event = std::make_shared<Reaction>(data);

    // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      // Check particle ID's and fill the reaction class
      if (data->id(part) == PIP) {
        event->SetPip(part);
      } else if (data->id(part) == PROTON) {
        event->SetProton(part);
      } else if (data->id(part) == PIM) {
        event->SetPim(part);
      } else {
        event->SetOther(part);
      }
    }

    // Check the reaction class what kind of even it is and fill the appropriate histograms
    if (event->TwoPion_missingPim()) {
      // _hists->WvsQ2_Fill(event->W(), event->Q2(), data->dc_sect(0));
      _hists->Fill_W_bin_check(event);
      _hists->Fill_histSevenD_pim(event);
      _hists->Fill_histSevenD_pip(event);
      _hists->Fill_histSevenD_prot(event);

   }
  }

  // Return the total number of events
  return total;
}

size_t run_files(std::vector<std::string> inputs, std::shared_ptr<Histogram> hists, int thread_id) {
  // Called once for each thread
  // Make a new chain to process for this thread
  auto chain = std::make_shared<TChain>("h10");
  // Add every file to the chain
  for (auto in : inputs) chain->Add(in.c_str());

  // Run the function over each thread
  return run(chain, hists, thread_id);
}

#endif
