//j
// Created by carlos on 13/05/19.
//

#ifndef BOOSTTOUR_ARC_H
#define BOOSTTOUR_ARC_H


class Arc {

private:
    int o, d, delay, jitter, bandwidth, estimateLinkDuration;

public:

    Arc(int o, int d, int delay, int jitter, int bandwidth, int estimateLinkDuration);

    int getO() const;

    void setO(int o);

    int getD() const;

    void setD(int d);

    int getDelay() const;

    void setDelay(int delay);

    int getJitter() const;

    void setJitter(int jitter);

    int getBandwidth() const;

    void setBandwidth(int bandwidth);

    int getEstimateLinkDuration() const;

    void setEstimateLinkDuration(int estimateLinkDuration);
};
#endif //BOOSTTOUR_ARC_H
