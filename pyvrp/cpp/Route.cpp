#include "Route.h"
#include "DurationSegment.h"
#include "LoadSegment.h"

#include <algorithm>
#include <fstream>
#include <numeric>

using pyvrp::Cost;
using pyvrp::Distance;
using pyvrp::Duration;
using pyvrp::Load;
using pyvrp::Route;

using Client = size_t;

Route::Route(ProblemData const &data, Visits visits, size_t const vehicleType)
    : visits_(std::move(visits)), centroid_({0, 0}), vehicleType_(vehicleType)
{
    auto const &vehType = data.vehicleType(vehicleType);
    startDepot_ = vehType.startDepot;
    endDepot_ = vehType.endDepot;

    DurationSegment ds = {startDepot_, vehType};
    std::vector<LoadSegment> loadSegments;
    loadSegments.reserve(data.numLoadDimensions());
    for (size_t i = 0; i != data.numLoadDimensions(); ++i)
        loadSegments.emplace_back(0, 0, 0);

    size_t prevClient = startDepot_;

    auto const &distances = data.distanceMatrix(vehType.profile);
    auto const &durations = data.durationMatrix(vehType.profile);

    for (size_t idx = 0; idx != size(); ++idx)
    {
        auto const client = visits_[idx];
        ProblemData::Client const &clientData = data.location(client);

        distance_ += distances(prevClient, client);
        travel_ += durations(prevClient, client);
        service_ += clientData.serviceDuration;
        prizes_ += clientData.prize;

        centroid_.first += static_cast<double>(clientData.x) / size();
        centroid_.second += static_cast<double>(clientData.y) / size();

        auto const clientDS = DurationSegment(client, clientData);
        ds = DurationSegment::merge(durations, ds, clientDS);

        for (size_t i = 0; i != data.numLoadDimensions(); ++i)
        {
            auto const clientLs = LoadSegment(clientData, i);
            loadSegments[i] = LoadSegment::merge(loadSegments[i], clientLs);
        }

        prevClient = client;
    }

    auto const last = visits_.empty() ? startDepot_ : visits_.back();
    distance_ += distances(last, endDepot_);
    distanceCost_ = vehType.unitDistanceCost * static_cast<Cost>(distance_);
    excessDistance_ = std::max<Distance>(distance_ - vehType.maxDistance, 0);

    travel_ += durations(last, endDepot_);

    delivery_.reserve(data.numLoadDimensions());
    pickup_.reserve(data.numLoadDimensions());
    excessLoad_.reserve(data.numLoadDimensions());
    for (size_t i = 0; i != data.numLoadDimensions(); ++i)
    {
        delivery_.push_back(loadSegments[i].delivery());
        pickup_.push_back(loadSegments[i].pickup());
        excessLoad_.push_back(
            std::max<Load>(loadSegments[i].load() - vehType.capacity[i], 0));
    }

    DurationSegment endDS(endDepot_, vehType);
    ds = DurationSegment::merge(durations, ds, endDS);
    duration_ = ds.duration();
    durationCost_ = vehType.unitDurationCost * static_cast<Cost>(duration_);
    startTime_ = ds.twEarly();
    slack_ = ds.twLate() - ds.twEarly();
    timeWarp_ = ds.timeWarp(vehType.maxDuration);
    release_ = ds.releaseTime();
}

Route::Route(Visits visits,
             Distance distance,
             Cost distanceCost,
             Distance excessDistance,
             std::vector<Load> delivery,
             std::vector<Load> pickup,
             std::vector<Load> excessLoad,
             Duration duration,
             Cost durationCost,
             Duration timeWarp,
             Duration travel,
             Duration service,
             Duration wait,
             Duration release,
             Duration startTime,
             Duration slack,
             Cost prizes,
             std::pair<double, double> centroid,
             size_t vehicleType,
             size_t startDepot,
             size_t endDepot)
    : visits_(std::move(visits)),
      distance_(distance),
      distanceCost_(distanceCost),
      excessDistance_(excessDistance),
      delivery_(std::move(delivery)),
      pickup_(std::move(pickup)),
      excessLoad_(std::move(excessLoad)),
      duration_(duration),
      durationCost_(durationCost),
      timeWarp_(timeWarp),
      travel_(travel),
      service_(service),
      wait_(wait),
      release_(release),
      startTime_(startTime),
      slack_(slack),
      prizes_(prizes),
      centroid_(centroid),
      vehicleType_(vehicleType),
      startDepot_(startDepot),
      endDepot_(endDepot)
{
}

size_t Route::numLoadDimensions() const { return excessLoad_.size(); }

bool Route::empty() const { return visits_.empty(); }

size_t Route::size() const { return visits_.size(); }

Client Route::operator[](size_t idx) const { return visits_[idx]; }

Route::Visits::const_iterator Route::begin() const { return visits_.cbegin(); }

Route::Visits::const_iterator Route::end() const { return visits_.cend(); }

Route::Visits const &Route::visits() const { return visits_; }

Distance Route::distance() const { return distance_; }

Cost Route::distanceCost() const { return distanceCost_; }

Distance Route::excessDistance() const { return excessDistance_; }

Load Route::delivery(size_t dimension) const
{
    if (dimension >= delivery_.size())
        throw std::out_of_range(
            "Dimension is out of range for the route's delivery load.");

    return delivery_[dimension];
}

Load Route::pickup(size_t dimension) const
{
    if (dimension >= pickup_.size())
        throw std::out_of_range(
            "Dimension is out of range for the route's pickup load.");

    return pickup_[dimension];
}

Load Route::excessLoad(size_t dimension) const
{
    if (dimension >= excessLoad_.size())
        throw std::out_of_range(
            "Dimension is out of range for the route's excess load.");

    return excessLoad_[dimension];
}

Duration Route::duration() const { return duration_; }

Cost Route::durationCost() const { return durationCost_; }

Duration Route::serviceDuration() const { return service_; }

Duration Route::timeWarp() const { return timeWarp_; }

Duration Route::waitDuration() const { return duration_ - travel_ - service_; }

Duration Route::travelDuration() const { return travel_; }

Duration Route::startTime() const { return startTime_; }

Duration Route::endTime() const { return startTime_ + duration_ - timeWarp_; }

Duration Route::slack() const { return slack_; }

Duration Route::releaseTime() const { return release_; }

Cost Route::prizes() const { return prizes_; }

std::pair<double, double> const &Route::centroid() const { return centroid_; }

size_t Route::vehicleType() const { return vehicleType_; }

size_t Route::startDepot() const { return startDepot_; }

size_t Route::endDepot() const { return endDepot_; }

bool Route::isFeasible() const
{
    return !hasExcessLoad() && !hasTimeWarp() && !hasExcessDistance();
}

bool Route::hasExcessLoad() const
{
    for (Load const &load : excessLoad_)
        if (load > 0)
            return true;

    return false;
}

bool Route::hasExcessDistance() const { return excessDistance_ > 0; }

bool Route::hasTimeWarp() const { return timeWarp_ > 0; }

bool Route::operator==(Route const &other) const
{
    // First compare simple attributes, since that's a quick and cheap check.
    // Only when these are the same we test if the visits are all equal.
    // clang-format off
    return distance_ == other.distance_
        && duration_ == other.duration_
        && timeWarp_ == other.timeWarp_
        && vehicleType_ == other.vehicleType_
        && visits_ == other.visits_;
    // clang-format on
}

std::ostream &operator<<(std::ostream &out, Route const &route)
{
    for (auto const client : route)
        out << client << ' ';
    return out;
}
