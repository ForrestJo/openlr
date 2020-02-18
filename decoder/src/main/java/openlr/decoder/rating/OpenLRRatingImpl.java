/**
 * Licensed to the TomTom International B.V. under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  TomTom International B.V.
 * licenses this file to you under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License.  You may obtain a copy of the License at
 * <p>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p>
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 * <p>
 * Copyright (C) 2009-2019 TomTom International B.V.
 * <p>
 * TomTom (Legal Department)
 * Email: legal@tomtom.com
 * <p>
 * TomTom (Technical contact)
 * Email: openlr@tomtom.com
 * <p>
 * Address: TomTom International B.V., Oosterdoksstraat 114, 1011DK Amsterdam,
 * the Netherlands
 * <p>
 */
package openlr.decoder.rating;

import openlr.LocationReferencePoint;
import openlr.OpenLRProcessingException;
import openlr.decoder.OpenLRDecoderProcessingException;
import openlr.decoder.properties.OpenLRDecoderProperties;
import openlr.map.*;
import openlr.map.utils.GeometryUtils;
import org.locationtech.jts.geom.Coordinate;

import java.util.Arrays;

import static openlr.decoder.OpenLRDecoderProcessingException.DecoderProcessingError.INVALID_MAP_DATA;

/**
 * The Class OpenLRRatingImpl implements the Frechet Distance function for OpenLR.
 * <p>
 * It is an experiment implementation to test how good the actual line value matches
 * the LRP attributes by Frechet Distance between them
 * <p>
 * OpenLR is a trade mark of TomTom International B.V.
 * <p>
 * email: software@openlr.org
 *
 * @author TomTom International B.V.
 */
public class OpenLRRatingImpl implements OpenLRRating {

    /** {@inheritDoc}
     * @return*/
    @Override
    public final int getRating(final OpenLRDecoderProperties properties,
                               final int distance, final LocationReferencePoint p,
                               final Line line, final int projectionAlongLine) throws OpenLRProcessingException {
        try {

            // Location reference point geometry
            GeoCoordinates lrpEndPoint = GeometryUtils.determineCoordinateInDistance(
                    p.getLongitudeDeg(),
                    p.getLatitudeDeg(),
                    (int) p.getBearing(),
                    (double) properties.getBearingDistance() / 1000);

            int offset = p.isLastLRP() ? -properties.getBearingDistance() : properties.getBearingDistance();

            // Points of the line to test geometry
            int endPointDistanceAlongLine = projectionAlongLine + offset;
            GeoCoordinates startPoint = line.getGeoCoordinateAlongLine(projectionAlongLine);
            GeoCoordinates endPoint = line.getGeoCoordinateAlongLine(endPointDistanceAlongLine);

            // Construct to coordinate array
            Coordinate lrpStartPointCoor = new Coordinate(p.getLongitudeDeg(), p.getLatitudeDeg());
            Coordinate lrpEndPointCoor = new Coordinate(lrpEndPoint.getLongitudeDeg(), lrpEndPoint.getLatitudeDeg());

            Coordinate startPointCoor = new Coordinate(startPoint.getLongitudeDeg(), startPoint.getLatitudeDeg());
            Coordinate endPointCoor = new Coordinate(endPoint.getLongitudeDeg(), endPoint.getLatitudeDeg());

            Coordinate[] lrpCoordinates = new Coordinate[]{lrpStartPointCoor, lrpEndPointCoor};
            Coordinate[] testLineCoordinates = new Coordinate[]{startPointCoor, endPointCoor};

            int lrpLength = lrpCoordinates.length;
            int testLength = testLineCoordinates.length;

            double[][] mem = new double[lrpLength][testLength];

            for (double[] row : mem)
                Arrays.fill(row, -1.0);

            double frechetDistance = computeDiscreteFrechetDistance(lrpCoordinates, testLineCoordinates, mem, lrpLength - 1, testLength - 1);

            double inverseFrechetDistance = (1 / frechetDistance);

            return (int) inverseFrechetDistance;
        }
        catch (InvalidMapDataException e) {
            throw new OpenLRDecoderProcessingException(INVALID_MAP_DATA, e);
        }
    }

    /**
     *
     * @param lrpCoordinates start and end coordinates of the Location reference point
     * @param testLineCoordinates start and end coordinates of the line to test
     * @param mem
     * @param i
     * @param j
     * @return Computed Discrete Frechet Distance in Double
     */

    private double computeDiscreteFrechetDistance(Coordinate[] lrpCoordinates, Coordinate[] testLineCoordinates, double[][] mem, int i, int j) {
        if (mem[i][j] > -1) {
            return mem[i][j];
        } else if (i == 0 && j == 0) {
            mem[i][j] = lrpCoordinates[i].distance(testLineCoordinates[j]);
        } else if (i > 0 && j == 0) {
            mem[i][j] = getMax(computeDiscreteFrechetDistance(lrpCoordinates, testLineCoordinates, mem, i - 1, j), lrpCoordinates[i].distance(testLineCoordinates[j]));
        } else if (i == 0 && j > 0) {
            mem[i][j] = getMax(computeDiscreteFrechetDistance(lrpCoordinates, testLineCoordinates, mem, i, j - 1), lrpCoordinates[i].distance(testLineCoordinates[j]));
        } else if (i > 0 && j > 0) {
            mem[i][j] = getMax(getMin(computeDiscreteFrechetDistance(lrpCoordinates, testLineCoordinates, mem, i - 1, j), computeDiscreteFrechetDistance(lrpCoordinates, testLineCoordinates, mem, i - 1, j - 1), computeDiscreteFrechetDistance(lrpCoordinates, testLineCoordinates, mem, i, j - 1)), lrpCoordinates[i].distance(testLineCoordinates[j]));
        } else {
            mem[i][j] = -2.0;
        }

        return mem[i][j];
    }

    private double getMax(double... values) {
        return Arrays.stream(values)
                .max()
                .getAsDouble();
    }

    private double getMin(double... values) {
        return Arrays.stream(values)
                .min()
                .getAsDouble();
    }
}
