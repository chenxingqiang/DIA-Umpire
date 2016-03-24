/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package MSUmpire.PeakDataStructure;

import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;

/**
 * B-spline smoothing
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class Bspline {

    private Bspline() {
    }

    private static class Tls {

        float[] arr_X = new float[1 << 5];
        float[] arr_Y = new float[1 << 5];
        float[] smoothed_X = new float[1 << 5];
        float[] smoothed_Y = new float[1 << 5];
    }

    private static final ThreadLocal<Tls> TLS = new ThreadLocal<Tls>() {
        @Override
        protected Tls initialValue() {
            return new Tls();
        }
    };

    private static Tls get_tls(final int n, final int num_pts) {
        final Tls t = TLS.get();
        if (t.arr_X.length < n) {
            t.arr_X = new float[n << 1];
            t.arr_Y = new float[n << 1];
        }
        if (t.smoothed_X.length < num_pts) {
            t.smoothed_X = new float[num_pts << 1];
            t.smoothed_Y = new float[num_pts << 1];
        }

        return t;
    }

    public static XYPointCollection Run__0(final XYPointCollection data, final int num_pts, final int smoothDegree) {
        assert smoothDegree==2;
        final int n = data.PointCount() + 2;
//    	final float[] arr_X=new float[n];
//    	final float[] arr_Y=new float[n];
//    	final float[] smoothed_X=new float[num_pts];
//    	final float[] smoothed_Y=new float[num_pts];
        final Tls t = get_tls(n, num_pts);

        final float[] arr_X = t.arr_X;
        final float[] arr_Y = t.arr_Y;
        final float[] smoothed_X = t.smoothed_X;
        final float[] smoothed_Y = t.smoothed_Y;
        {
            final XYData pt = data.Data.get(0);
            arr_X[0] = arr_X[1] = pt.getX();
            arr_Y[0] = arr_Y[1] = pt.getY();
        }
        for (int i = 1; i < n - 3; ++i) {
            final XYData pt = data.Data.get(i);
            arr_X[i + 1] = pt.getX();
            arr_Y[i + 1] = pt.getY();
        }
        {
            final XYData pt = data.Data.get(n - 3);
            arr_X[n - 2] = arr_X[n - 1] = pt.getX();
            arr_Y[n - 2] = arr_Y[n - 1] = pt.getY();
        }
        QuadBspline.get_smoothed(arr_X, n, num_pts, smoothed_X);
        QuadBspline.get_smoothed(arr_Y, n, num_pts, smoothed_Y);
        final XYPointCollection smoothed = new XYPointCollection();
        for (int i = 0; i < num_pts; ++i) {
            smoothed.AddPoint(smoothed_X[i], smoothed_Y[i]);
        }
        smoothed.Data.Finalize();
        return smoothed;
    }
    public static XYPointCollection Run(final XYPointCollection data, int num_pts, final int smoothDegree) {
        // this follows closely CC's implementation of bspline
        assert smoothDegree==2;
        ++num_pts; // to follow CC's method
        final int n = data.PointCount();
        final Tls t = get_tls(n, num_pts);

        final float[] arr_X = t.arr_X;
        final float[] arr_Y = t.arr_Y;
        final float[] smoothed_X = t.smoothed_X;
        final float[] smoothed_Y = t.smoothed_Y;

        for (int i = 0; i < n; ++i) {
            final XYData pt = data.Data.get(i);
            arr_X[i] = pt.getX();
            arr_Y[i] = pt.getY();
        }
        QuadBspline.get_smoothed(arr_X, n, num_pts, smoothed_X);
        QuadBspline.get_smoothed(arr_Y, n, num_pts, smoothed_Y);
        final XYPointCollection smoothed = new XYPointCollection();
        for (int i = 0; i < num_pts; ++i) {
            smoothed.AddPoint(smoothed_X[i], smoothed_Y[i]);
        }
        smoothed.Data.Finalize();
        return smoothed;
    }

}

class QuadBspline {
    //https://en.wikipedia.org/wiki/B-spline

    private QuadBspline() {
        throw new AssertionError();
    }

    private static float square(final float num) {
        return num * num;
    }

    /**
     * first two entries of data must be equal last two entries of data must be
     * equal
     *
     * @param data
     * @param num_pts
     * @return
     */
    public static void get_smoothed(final float[] data, final int n,
            final int num_pts, final float[] ret) {

//        assert data[0] == data[1];
//        assert data[n - 1] == data[n - 2];
        final int p = 2;//quadratic spline
        final float int_len = ((n + p - 2) - 2) / (float) (num_pts - 1);

        ret[0] = data[0];
        for (int i = 1; i < num_pts - 1; ++i) {
            final float t = 2 + i * int_len;
            final int ti = (int) t;
            final float rem = t - ti;
            final float b0 = square(1 - rem) / 2,
                    b2 = square(rem) / 2,
                    b1 = 1 - b0 - b2;
            ret[i] = b0 * data[ti - 2] + b1 * data[ti - 1] + b2 * data[ti];
        }
        ret[num_pts - 1] = data[n - 1];
    }
}
