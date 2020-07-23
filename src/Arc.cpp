//
// Created by carlos on 13/05/19.
//

#include "Arc.h"


Arc::Arc(int o, int d, int delay, int jitter, int bandwidth, int estimateLinkDuration) : o(o), d(d), delay(delay),
                                                                                         jitter(jitter),
                                                                                         bandwidth(bandwidth),
                                                                                         estimateLinkDuration(
                                                                                                 estimateLinkDuration) {}

int Arc::getO() const {
    return o;
}

void Arc::setO(int o) {
    Arc::o = o;
}

int Arc::getD() const {
    return d;
}

void Arc::setD(int d) {
    Arc::d = d;
}

int Arc::getDelay() const {
    return delay;
}

void Arc::setDelay(int delay) {
    Arc::delay = delay;
}

int Arc::getJitter() const {
    return jitter;
}

void Arc::setJitter(int jitter) {
    Arc::jitter = jitter;
}

int Arc::getBandwidth() const {
    return bandwidth;
}

void Arc::setBandwidth(int bandwidth) {
    Arc::bandwidth = bandwidth;
}

int Arc::getEstimateLinkDuration() const {
    return estimateLinkDuration;
}

void Arc::setEstimateLinkDuration(int estimateLinkDuration) {
    Arc::estimateLinkDuration = estimateLinkDuration;
}

