/** $lic$
 * Copyright (C) 2012-2015 by Massachusetts Institute of Technology
 * Copyright (C) 2010-2013 by The Board of Trustees of Stanford University
 *
 * This file is part of zsim.
 *
 * zsim is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, version 2.
 *
 * If you use this software in your research, we request that you reference
 * the zsim paper ("ZSim: Fast and Accurate Microarchitectural Simulation of
 * Thousand-Core Systems", Sanchez and Kozyrakis, ISCA-40, June 2013) as the
 * source of the simulator in any publications that use this software, and that
 * you send us a citation of your work.
 *
 * zsim is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "coherence_ctrls.h"
#include "cache.h"
#include "network.h"

#include "zsim.h"
#include "event_recorder.h"
#include "timing_event.h"

/* Do a simple XOR block hash on address to determine its bank. Hacky for now,
 * should probably have a class that deals with this with a real hash function
 * (TODO)
 */
uint32_t MESIBottomCC::getParentId(Address lineAddr) {
    //Hash things a bit
    uint32_t res = 0;
    uint64_t tmp = lineAddr;
    for (uint32_t i = 0; i < 4; i++) {
        res ^= (uint32_t) ( ((uint64_t)0xffff) & tmp);
        tmp = tmp >> 16;
    }
    return (res % parents.size());
}


void MESIBottomCC::init(const g_vector<MemObject*>& _parents, Network* network, const char* name) {
    parents.resize(_parents.size());
    parentRTTs.resize(_parents.size());
    for (uint32_t p = 0; p < parents.size(); p++) {
        parents[p] = _parents[p];
        parentRTTs[p] = (network)? network->getRTT(name, parents[p]->getName()) : 0;
    }
}


uint64_t MESIBottomCC::processEviction(Address wbLineAddr, uint32_t lineId, bool lowerLevelWriteback, uint64_t cycle, uint32_t srcId) {
    MESIState* state = &array[lineId];
    if (lowerLevelWriteback) {
        //If this happens, when tcc issued the invalidations, it got a writeback. This means we have to do a PUTX, i.e. we have to transition to M if we are in E
        assert(*state == M || *state == E); //Must have exclusive permission!
        *state = M; //Silent E->M transition (at eviction); now we'll do a PUTX
    }
    uint64_t respCycle = cycle;
    switch (*state) {
        case I:
            break; //Nothing to do
        case S:
        case E:
            {
                MemReq req = {wbLineAddr, PUTS, selfId, state, cycle, &ccLock, *state, srcId, 0 /*no flags*/};
                respCycle = parents[getParentId(wbLineAddr)]->access(req);
            }
            break;
        case M:
            {
                MemReq req = {wbLineAddr, PUTX, selfId, state, cycle, &ccLock, *state, srcId, 0 /*no flags*/};
                respCycle = parents[getParentId(wbLineAddr)]->access(req);
            }
            break;

        default: panic("!?");
    }
    assert_msg(*state == I, "Wrong final state %s on eviction", MESIStateName(*state));
    return respCycle;
}

uint64_t MESIBottomCC::processAccess(Address lineAddr, uint32_t lineId, AccessType type, uint64_t cycle, uint32_t srcId, uint32_t flags) {
    uint64_t respCycle = cycle;
    MESIState* state = &array[lineId];
    switch (type) {
        // A PUTS/PUTX does nothing w.r.t. higher coherence levels --- it dies here
        case PUTS: //Clean writeback, nothing to do (except profiling)
            assert(*state != I);
            profPUTS.inc();
            break;
        case PUTX: //Dirty writeback
            assert(*state == M || *state == E);
            if (*state == E) {
                //Silent transition, record that block was written to
                *state = M;
            }
            profPUTX.inc();
            break;
        case GETS:
            if (*state == I) {
                uint32_t parentId = getParentId(lineAddr);
                MemReq req = {lineAddr, GETS, selfId, state, cycle, &ccLock, *state, srcId, flags};
                uint32_t nextLevelLat = parents[parentId]->access(req) - cycle;
                uint32_t netLat = parentRTTs[parentId];
                profGETNextLevelLat.inc(nextLevelLat);
                profGETNetLat.inc(netLat);
                respCycle += nextLevelLat + netLat;
                profGETSMiss.inc();
                assert(*state == S || *state == E);
            } else {
                profGETSHit.inc();
            }
            break;
        case GETX:
            if (*state == I || *state == S) {
                //Profile before access, state changes
                if (*state == I) profGETXMissIM.inc();
                else profGETXMissSM.inc();
                uint32_t parentId = getParentId(lineAddr);
                MemReq req = {lineAddr, GETX, selfId, state, cycle, &ccLock, *state, srcId, flags};
                uint32_t nextLevelLat = parents[parentId]->access(req) - cycle;
                uint32_t netLat = parentRTTs[parentId];
                profGETNextLevelLat.inc(nextLevelLat);
                profGETNetLat.inc(netLat);
                respCycle += nextLevelLat + netLat;
            } else {
                if (*state == E) {
                    // Silent transition
                    // NOTE: When do we silent-transition E->M on an ML hierarchy... on a GETX, or on a PUTX?
                    /* Actually, on both: on a GETX b/c line's going to be modified anyway, and must do it if it is the L1 (it's OK not
                     * to transition if L2+, we'll TX on the PUTX or invalidate, but doing it this way minimizes the differences between
                     * L1 and L2+ controllers); and on a PUTX, because receiving a PUTX while we're in E indicates the child did a silent
                     * transition and now that it is evictiong, it's our turn to maintain M info.
                     */
                    *state = M;
                }
                profGETXHit.inc();
            }
            assert_msg(*state == M, "Wrong final state on GETX, lineId %d numLines %d, finalState %s", lineId, numLines, MESIStateName(*state));
            break;

        default: panic("!?");
    }
    assert_msg(respCycle >= cycle, "XXX %ld %ld", respCycle, cycle);
    return respCycle;
}

void MESIBottomCC::processWritebackOnAccess(Address lineAddr, uint32_t lineId, AccessType type) {
    MESIState* state = &array[lineId];
    assert(*state == M || *state == E);
    if (*state == E) {
        //Silent transition to M if in E
        *state = M;
    }
}

void MESIBottomCC::processInval(Address lineAddr, uint32_t lineId, InvType type, bool* reqWriteback) {
    MESIState* state = &array[lineId];
    assert(*state != I);
    switch (type) {
        case INVX: //lose exclusivity
            //Hmmm, do we have to propagate loss of exclusivity down the tree? (nah, topcc will do this automatically -- it knows the final state, always!)
            assert_msg(*state == E || *state == M, "Invalid state %s", MESIStateName(*state));
            if (*state == M) *reqWriteback = true;
            *state = S;
            profINVX.inc();
            break;
        case INV: //invalidate
            assert(*state != I);
            if (*state == M) *reqWriteback = true;
            *state = I;
            profINV.inc();
            break;
        case FWD: //forward
            assert_msg(*state == S, "Invalid state %s on FWD", MESIStateName(*state));
            profFWD.inc();
            break;
        default: panic("!?");
    }
    //NOTE: BottomCC never calls up on an invalidate, so it adds no extra latency
}


uint64_t MESIBottomCC::processNonInclusiveWriteback(Address lineAddr, AccessType type, uint64_t cycle, MESIState* state, uint32_t srcId, uint32_t flags) {
    if (!nonInclusiveHack) panic("Non-inclusive %s on line 0x%lx, this cache should be inclusive", AccessTypeName(type), lineAddr);

    //info("Non-inclusive wback, forwarding");
    MemReq req = {lineAddr, type, selfId, state, cycle, &ccLock, *state, srcId, flags | MemReq::NONINCLWB};
    uint64_t respCycle = parents[getParentId(lineAddr)]->access(req);
    return respCycle;
}


/* MESITopCC implementation */

void MESITopCC::init(const g_vector<BaseCache*>& _children, Network* network, const char* name, MESICC* cc) {
    if (_children.size() > MAX_CACHE_CHILDREN) {
        panic("[%s] Children size (%d) > MAX_CACHE_CHILDREN (%d)", name, (uint32_t)_children.size(), MAX_CACHE_CHILDREN);
    }
    children.resize(_children.size());
    childrenRTTs.resize(_children.size());
    for (uint32_t c = 0; c < children.size(); c++) {
        children[c] = _children[c];
        childrenRTTs[c] = (network)? network->getRTT(name, children[c]->getName()) : 0;
    }
    this->cc = cc;
}

uint64_t MESITopCC::sendInvalidates(Address lineAddr, uint32_t lineId, InvType type, bool* reqWriteback, uint64_t cycle, uint32_t srcId) {
    //Send down downgrades/invalidates
    Entry* e = &array[lineId];

    //Don't propagate downgrades if sharers are not exclusive.
    if (type == INVX && !e->isExclusive()) {
        return cycle;
    }

    uint64_t maxCycle = cycle; //keep maximum cycle only, we assume all invals are sent in parallel
    if (!e->isEmpty()) {
        uint32_t numChildren = children.size();
        uint32_t sentInvs = 0;
        for (uint32_t c = 0; c < numChildren; c++) {
            if (e->sharers[c]) {
                InvReq req = {lineAddr, type, reqWriteback, cycle, srcId};
                uint64_t respCycle = children[c]->invalidate(req);
                respCycle += childrenRTTs[c];
                maxCycle = MAX(respCycle, maxCycle);
                if (type == INV) e->sharers[c] = false;
                sentInvs++;
            }
        }
        assert(sentInvs == e->numSharers);
        if (type == INV) {
            e->numSharers = 0;
        } else {
            //TODO: This is kludgy -- once the sharers format is more sophisticated, handle downgrades with a different codepath
            assert(e->exclusive);
            assert(e->numSharers == 1);
            e->exclusive = false;
        }
    }
    return maxCycle;
}


uint64_t MESITopCC::processEviction(Address wbLineAddr, uint32_t lineId, bool* reqWriteback, uint64_t cycle, uint32_t srcId) {
    if (nonInclusiveHack) {
        // Don't invalidate anything, just clear our entry
        array[lineId].clear();
        return cycle;
    } else {
        //Send down invalidates
        return sendInvalidates(wbLineAddr, lineId, INV, reqWriteback, cycle, srcId);
    }
}

uint64_t MESITopCC::processAccess(Address lineAddr, uint32_t lineId, AccessType type, uint32_t childId, bool haveExclusive,
                                  MESIState* childState, bool* inducedWriteback, uint64_t cycle, uint32_t srcId, uint32_t flags, CacheArray* data_array, MemReq req) {
    Entry* e = &array[lineId];
    
    std::vector<uint32_t> id_vector;
    std::vector<Address> address_vector;
    uint32_t assoc = data_array->get_set(lineAddr, id_vector, address_vector);
    //info("data_array->get_set() = %d\n", assoc);
    //if (e->numSharers) info("e->numSharers = %d\n", e->numSharers);

    assert(id_vector.size() == assoc);
    int32_t way = -1; // stores the way that the current cache line is in=
    for(uint32_t i = 0; i < id_vector.size(); i++){
        if (id_vector[i] == lineId){
            way = i;
        }
    }
    assert(way != -1);

    // Ways 20 and 21 (the last two) are the reclaimed cache lines)
    bool accessing_reclaimed_cache_line = false;
    if (enable_reclaim && ((way == 20) || (way == 21))){
        accessing_reclaimed_cache_line = true;
    }

    uint32_t min_way = (way/10)*10; // stores the starting way that we need to check (0, 10, 20)
    way = way + 1;
    uint32_t max_way = (way/10 + (way % 10 != 0))*10 - 1; // stores the last way that we need to check (9, 19, 29)

    // This variable says if we are one of the cache lines whose directories can be repurposed
    bool reclaimable_way = false;
    if ((max_way < assoc) && (assoc > 8)){ // assoc > 8 means we are the LLC - won't work for all configs
        reclaimable_way = true;
    }

    bool reclaiming_possible_before = true;
    if (reclaimable_way){ 
        for(uint32_t i = min_way; i <= max_way; i++){
            if (array[id_vector[i]].numSharers > 1){
                reclaiming_possible_before = false;
            }
        }
    }

    //info("reclaiming_possible = %d", reclaiming_possible_before);

    uint64_t respCycle = cycle;
    switch (type) {
        case PUTX:
            assert(e->isExclusive());
            if (flags & MemReq::PUTX_KEEPEXCL) {
                assert(e->sharers[childId]);
                assert(*childState == M);
                *childState = E; //they don't hold dirty data anymore
                break; //don't remove from sharer set. It'll keep exclusive perms.
            }
            //note NO break in general
        case PUTS:
            assert(e->sharers[childId]);
            e->sharers[childId] = false;
            e->numSharers--;
            *childState = I;
            break;
        case GETS:
            if (e->isEmpty() && haveExclusive && !(flags & MemReq::NOEXCL)) {
                //Give in E state
                e->exclusive = true;
                e->sharers[childId] = true;
                e->numSharers = 1;
                *childState = E;
            } else {
                //Give in S state
                assert(e->sharers[childId] == false);

                bool need_lim_ptr_eviction = false;
                if (e->isExclusive()) {
                    //Downgrade the exclusive sharer
                    if (accessing_reclaimed_cache_line){
                        // if we are a reclaimed cache line, we need to completely invalidate the exclusive copy
                        // as we can only hold a single sharer
                        respCycle = sendInvalidates(lineAddr, lineId, INV, inducedWriteback, cycle, srcId);
                    }
                    else{
                        respCycle = sendInvalidates(lineAddr, lineId, INVX, inducedWriteback, cycle, srcId);
                    }
                    //respCycle = sendInvalidates(lineAddr, lineId, INVX, inducedWriteback, cycle, srcId);
                }
                else{
                    // CS533 
                    // LIMITED POINTER IMPLEMENTATION
                    uint32_t MAX_SHARERS;
                    if (accessing_reclaimed_cache_line){
                        MAX_SHARERS = 1;
                    }
                    else{
                        MAX_SHARERS = max_sharers;
                    }
                    
                    if (MAX_SHARERS != 0){
                        // This code functionally implements the lim. pointer scheme
                        // while still storing sharer info in the full bit vector
                        // that zsim uses
                        assert(e->numSharers <= MAX_SHARERS);
                        if (e->numSharers == MAX_SHARERS){
                            need_lim_ptr_eviction = true;
                            assert(e->exclusive == false);

                            // send eviction to a random sharer
                            uint32_t rand_sharer = rand() % MAX_SHARERS;
                            uint32_t numChildren = children.size();
                            uint32_t sharer_ctr = 0;
                            bool sent_invalidation = false;
                            for (uint32_t c = 0; c < numChildren; c++){
                                if (e->sharers[c]){
                                    if (sharer_ctr == rand_sharer){
                                        //bool false_var = false;
                                        //InvReq req = {lineAddr, INV, &false_var, cycle, srcId};
                                        InvReq req = {lineAddr, INV, inducedWriteback, cycle, srcId};
                                        respCycle = children[c]->invalidate(req);
                                        respCycle += childrenRTTs[c];
                                        respCycle = MAX(respCycle, cycle);
                                        e->sharers[c] = false;
                                        sent_invalidation = true;
                                        break;
                                    }
                                    
                                    sharer_ctr++;
                                }        
                            }

                            assert(sent_invalidation);
                        }
                    }
                }

                assert_msg(!e->isExclusive(), "Can't have exclusivity here. isExcl=%d excl=%d numSharers=%d", e->isExclusive(), e->exclusive, e->numSharers);

                if (!need_lim_ptr_eviction) e->numSharers++;

                e->sharers[childId] = true;
                
                e->exclusive = false; //dsm: Must set, we're explicitly non-exclusive
                *childState = S;
            }
            break;
        case GETX:
            assert(haveExclusive); //the current cache better have exclusive access to this line

            // If child is in sharers list (this is an upgrade miss), take it out
            if (e->sharers[childId]) {
                assert_msg(!e->isExclusive(), "Spurious GETX, childId=%d numSharers=%d isExcl=%d excl=%d", childId, e->numSharers, e->isExclusive(), e->exclusive);
                e->sharers[childId] = false;
                e->numSharers--;
            }

            // Invalidate all other copies
            respCycle = sendInvalidates(lineAddr, lineId, INV, inducedWriteback, cycle, srcId);

            // Set current sharer, mark exclusive
            e->sharers[childId] = true;
            e->numSharers++;
            e->exclusive = true;

            assert(e->numSharers == 1);

            *childState = M; //give in M directly
            break;

        default: panic("!?");
    }

    bool reclaiming_possible_after = true;
    if (reclaimable_way){
        for(uint32_t i = min_way; i <= max_way; i++){
            if (array[id_vector[i]].numSharers > 1){
                reclaiming_possible_after = false;
            }
        }

        if (reclaiming_possible_before == false && reclaiming_possible_after == true){
            // increment number of reclaimable lines
            //info("reclaiming N -> Y");
            profNonReclaimableLines.dec();
            //info("num_reclaimed_sets = %ld", num_reclaimed_sets);
        }
        if (reclaiming_possible_before == true && reclaiming_possible_after == false){
            // decrement number of reclaimable lines
            //info("reclaiming Y -> N");
            profNonReclaimableLines.inc();

            // Enforce single-record invariant: Writeback access may have a timing
            // record. If so, read it.

            // Trying to fix a bug dealing with the timing records.
            // Usually, processAccess can only generate on timing record, which basically happens on a GETS or GETX, and there is a cache miss
            
            // A timing record can also be generated from the eviction, when we miss in the cache and allocate a new line

            // I don't care about very precisely modelling contention, etc. but it seems that the timing stuff is important for zsim

            /*
            EventRecorder* evRec = zinfo->eventRecorders[req.srcId];
            TimingRecord wbAcc;
            wbAcc.clear();
            if (unlikely(evRec && evRec->hasRecord())) {
                //info("TIMING RECORD!");

                wbAcc = evRec->popRecord();
            }*/

            //*
            EventRecorder* evRec = zinfo->eventRecorders[req.srcId];
            TimingRecord rec;
            rec.clear();
            bool hadRec = false;
            if (evRec && evRec->hasRecord()){
                hadRec = true;
                rec = evRec->popRecord();
            }//*/

            if (max_way == 9){
                //req.srcId = id_vector[20];
                // we don't have access to the access latency in this module, and I'm too lazy to add it
                // so i'm just adding 40 to respCycle (lower numbers don't work..........)
                if (cc->bcc->array[id_vector[20]] != I){
                    info("DOING EVICTION OF RECLAIMED LINE 1");
                    info("id_vector[20] == %d", id_vector[20]);
                    info("address_vector[20] == %lx", address_vector[20]);
                    cc->processEviction(req, address_vector[20], id_vector[20], respCycle+60);
                }
            }
            else if (max_way == 19){
                //req.srcId = id_vector[21];
                //cc->processEviction(req, address_vector[21], id_vector[21], respCycle+30);
                
                if (cc->bcc->array[id_vector[21]] != I){
                    info("DOING EVICTION OF RECLAIMED LINE 2");
                    info("id_vector[21] == %d", id_vector[21]);
                    info("address_vector[21] == %lx", address_vector[21]);
                    cc->processEviction(req, address_vector[21], id_vector[21], respCycle+60);
                }
            }

            if (evRec && evRec->hasRecord()){
                evRec->popRecord();
            }
            if (hadRec){
                evRec->pushRecord(rec);
            }
            
            // /*/
            //respCycle = cc->processAccess(req, lineId, respCycle, array);
            /*
            if (unlikely(wbAcc.isValid())) {
                if (!evRec->hasRecord()) {
                    //info("(!evRec->hasRecord())");
                    // Downstream should not care about endEvent for PUTs
                    wbAcc.endEvent = nullptr;
                    evRec->pushRecord(wbAcc);
                } else {
                    
                    //info("else");
                    //info("req.cycle = %ld", req.cycle);
                    // Connect both events
                    TimingRecord acc = evRec->popRecord();
                    //info("wbAcc.reqCycle = %ld", wbAcc.reqCycle);
                    assert(wbAcc.reqCycle >= req.cycle);
                    //info("acc.reqCycle = %ld", acc.reqCycle);
                    assert(acc.reqCycle >= req.cycle);
                    DelayEvent* startEv = new (evRec) DelayEvent(0);
                    DelayEvent* dWbEv = new (evRec) DelayEvent(wbAcc.reqCycle - req.cycle);
                    DelayEvent* dAccEv = new (evRec) DelayEvent(acc.reqCycle - req.cycle);

                    //info("acc.reqCycle = %ld", acc.reqCycle);
                    startEv->setMinStartCycle(req.cycle);
                    dWbEv->setMinStartCycle(req.cycle);
                    dAccEv->setMinStartCycle(req.cycle);
                    startEv->addChild(dWbEv, evRec)->addChild(wbAcc.startEvent, evRec);
                    startEv->addChild(dAccEv, evRec)->addChild(acc.startEvent, evRec);

                    acc.reqCycle = req.cycle;
                    acc.startEvent = startEv;
                    // endEvent / endCycle stay the same; wbAcc's endEvent not connected
                    evRec->pushRecord(acc);
                }
            }*/
        }
    }
    
    return respCycle;
}

uint64_t MESITopCC::processInval(Address lineAddr, uint32_t lineId, InvType type, bool* reqWriteback, uint64_t cycle, uint32_t srcId) {
    if (type == FWD) {//if it's a FWD, we should be inclusive for now, so we must have the line, just invLat works
        assert(!nonInclusiveHack); //dsm: ask me if you see this failing and don't know why
        return cycle;
    } else {
        //Just invalidate or downgrade down to children as needed
        return sendInvalidates(lineAddr, lineId, type, reqWriteback, cycle, srcId);
    }
}

