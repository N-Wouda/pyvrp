#ifndef HGS_PROBLEMDATA_H
#define HGS_PROBLEMDATA_H

#include "Matrix.h"
#include "Measure.h"
#include "XorShift128.h"

#include <iosfwd>
#include <vector>

class ProblemData
{
public:
    struct Client
    {
        int const x;                     // Coordinate X
        int const y;                     // Coordinate Y
        Capacity const demand;           // Demand
        Duration const serviceDuration;  // Service duration
        Duration const twEarly;          // Earliest possible arrival
        Duration const twLate;           // Latest possible arrival
        Cost const prize = 0;            // Prize for visiting this client
        bool const required = true;      // Must client be in solution?

        Client(int x,
               int y,
               int demand = 0,
               int serviceDuration = 0,
               int twEarly = 0,
               int twLate = 0,
               int prize = 0,
               bool required = true);
    };

private:
    Matrix<Distance> const dist_;  // Distance matrix (+depot)
    Matrix<Duration> const dur_;   // Duration matrix (+depot)
    std::vector<Client> clients_;  // Client (+depot) information

    size_t const numClients_;
    size_t const numVehicles_;
    Capacity const vehicleCapacity_;

public:
    /**
     * @param client Client whose data to return.
     * @return A struct containing the indicated client's information.
     */
    [[nodiscard]] inline Client const &client(size_t client) const;

    /**
     * @return A struct containing the depot's information.
     */
    [[nodiscard]] Client const &depot() const;

    /**
     * Returns the distance between the indicated two clients.
     *
     * @param first  First client.
     * @param second Second client.
     * @return distance from the first to the second client.
     */
    [[nodiscard]] inline Distance dist(size_t first, size_t second) const;

    /**
     * Returns the travel duration between the indicated two clients.
     *
     * @param first  First client.
     * @param second Second client.
     * @return Travel duration from the first to the second client.
     */
    [[nodiscard]] inline Duration duration(size_t first, size_t second) const;

    /**
     * @return The full travel distance matrix.
     */
    [[nodiscard]] Matrix<Distance> const &distanceMatrix() const;

    /**
     * @return The full travel duration matrix.
     */
    [[nodiscard]] Matrix<Duration> const &durationMatrix() const;

    /**
     * @return Total number of clients in this instance.
     */
    [[nodiscard]] size_t numClients() const;

    /**
     * @return Total number of vehicles available in this instance.
     */
    [[nodiscard]] size_t numVehicles() const;

    /**
     * @return Capacity of each vehicle in this instance.
     */
    [[nodiscard]] Capacity vehicleCapacity() const;

    /**
     * Constructs a ProblemData object with the given data. Assumes the list of
     * clients contains the depot, such that each vector is one longer than the
     * number of clients.
     *
     * @param clients      List of clients (including depot at index 0).
     * @param numVehicles  Number of vehicles.
     * @param vehicleCap   Vehicle capacity.
     * @param distMat      Distance matrix.
     * @param durMat       Duration matrix.
     */
    ProblemData(std::vector<Client> const &clients,
                size_t numVehicles,
                Capacity vehicleCap,
                Matrix<Distance> const distMat,
                Matrix<Duration> const durMat);
};

ProblemData::Client const &ProblemData::client(size_t client) const
{
    return clients_[client];
}

Distance ProblemData::dist(size_t first, size_t second) const
{
    return dist_(first, second);
}

Duration ProblemData::duration(size_t first, size_t second) const
{
    return dur_(first, second);
}

#endif  // HGS_PROBLEMDATA_H
