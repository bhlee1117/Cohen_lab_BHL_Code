/*
* arrayProduct.cpp - example in MATLAB External Interfaces
*
* Multiplies an input scalar (multiplier)
* times a MxN matrix (inMatrix)
* and outputs a MxN matrix (outMatrix)
*
* Usage : from MATLAB
*         >> outMatrix = arrayProduct(multiplier, inMatrix)
*
* This is a C++ MEX-file for MATLAB.
* Copyright 2017 The MathWorks, Inc.
*
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <cstdlib>
#include <cmath>
#include <complex>

class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        checkArguments(outputs, inputs);
        matlab::data::TypedArray<std::complex<double>> sptField = std::move(inputs[0]);
        matlab::data::TypedArray<double> sptPhOffset = std::move(inputs[1]);
        matlab::data::TypedArray<std::complex<double>> input_intensity = std::move(inputs[2]);
        matlab::data::TypedArray<double> wt = std::move(inputs[3]);
        double iterations = inputs[4][0];
        
        
        
        slmPh = wgs_loop(sptField,sptPhOffset,input_intensity,wt,iterations );
        
        
        outputs[0] = std::move(slmPh);
    }
/*------------------- wgs_loop -------------------*/
    matlab::data::TypedArray<double> wgs_loop(
            matlab::data::TypedArray<std::complex<double>>& sptField, 
            matlab::data::TypedArray<double>& sptPhOffset, 
            matlab::data::TypedArray<std::complex<double>>& input_intensity, 
            matlab::data::TypedArray<double>& wt, 
            double iterations) {
       
        for (i=0; i<30; i++) {
//             slmPh = angle(slmField);
//             sptField = sum(input_intensity.*exp(1i*(-slmPh+sptPh)),[1 2])/numel(xj);
//             wt = mean(abs(sptField))/abs(sptField).*wt;
//             sptPhOffset = angle(sptField);%./sptField;
//             slmField = sum(exp(1i*(-sptPhOffset-sptPh)).*wt,3);
        }
        
        return angle(slmField);
    }
/*------------------- wgs_loop END -------------------*/
    

/*------------------- add_diff_size -------------------*/
    matlab::data::TypedArray<double> add_diff_size(
            matlab::data:TypedArray<std::complex<double>> in1, // 2D array
            matlab::data:TypedArray<std::complex<double>> in2) // 3D array
    {
        using namespace matlab::data;
        ArrayFactory factory;
        
        std::vector<size_t,3> mat_dim{
                                    in1.getDimensions()[0],
                                    in1.getDimensions()[1],
                                    in2.getDimensions()[2]}
//         mat_dim[0] = in1.getDimensions()[0];
//         mat_dim[1] = in1.getDimensions()[1];
//         mat_dim[2] = in2.getDimensions()[2];
        if (in1.getDimensions().sizeof() == in2.getDimensions().sizeof())
        {
            auto out = factory.createArray<std::complex<double>>(in1.getDimensions());
            for (int i=0; i < in1.getNumElements()[0]; i++)
            {
                out[i] = in1[i]+in2[i];
            }
        }
        
        if (in2.getDimensions()[0] < in1.getDimensions()[0])
        {
            auto in2 = factory.createArray<std::complex<double>>(mat_dim);
            for (i = 0; i < mat_dim[0]; i++)
            {
                for (j = 0; j < mat_dim[1]; j++)
                {
                    for (k = 0; k < mat_dim[2]; k++)
                    {
                        in2_temp[i][j][k] = in1[i][j];
                    }
                }
            }
        }
        auto in1 = factory.createArray<std::complex<double>>(mat_dim);
        for (i = 0; i < mat_dim[0]; i++)
        {
            for (j = 0; j < mat_dim[1]; j++)
            {
                for (k = 0; k < mat_dim[2]; k++)
                {
                    in1[i][j][k] = in1[i][j];
                }
            }
        }
        
        auto out = factory.createArray<std::complex<double>>(mat_dim);
        for (i = 0; i < mat_dim[0]; i++)
        {
            for (j = 0; j < mat_dim[1]; j++)
            {
                for (k = 0; k < mat_dim[2]; k++)
                {
                    out[i][j][k] = in1[i][j][k] + in2[i][j][k];
                }
            }
        }
        
        return out;
    }
/*------------------- add_diff_size END-------------------*/
            
    
/*------------------- divide_diff_size -------------------*/
    
    
            
/*------------------- divide_diff_size END -------------------*/
            
    
/*------------------- exp_array -------------------*/
    matlab::data::TypedArray<std::complex<double>> exp_array(matlab::data:TypedArray<std::complex<double>> in) { 
        using namespace matlab::data;
        ArrayFactory factory;
        
        ArrayDimension in_size = in.getDimensions();
        
        auto out = factory.createArray<std::complex<double>>(in.getDimensions);
                
        for (i_rows = 0; i_rows < in_size[0]; i_rows++){
            for (i_cols = 0; i_cols < in_size[1]; i_cols++){
                out[i_rows][i_cols] = std::exp(1i * in[i_rows][i_cols]);
            }
        }
        
        return out;
    }
/*------------------- exp_array END -------------------*/
    
    
/*------------------- arg_array -------------------*/
    matlab::data::TypedArray<double> arg_array(matlab::data:TypedArray<std::complex<double>> in)
    {
        using namespace matlab::data;
        ArrayFactory factory;
        TypedArray<double> out = factory.createArray<double>(in.getDimensions);
        
        ArrayDimension in_size = in.getDimensions();
        
        for (i_rows = 0; i_rows < in_size[0]; i_rows++){
            for (i_cols = 0; i_cols < in_size[1]; i_cols++){
                out[i_rows][i_cols] = std::arg(in[i_rows][i_cols]);
            }
        }
        
        return out;
    }
/*------------------- arg_array END -------------------*/
    
   
            
/*------------------- sum_array_01 -------------------*/
    auto sum_array_01(matlab::data::Array in){
        using namespace matlab::data;
        ArrayFactory factory;
        auto out = factory.createArray<std::complex<double>>(in.getDimensions()[2]);
        auto idx_lim = in.getDimensions();
        for (k = 0; k < idx_lim[2]; k++)
        {
            out[k] = 0;
            for (i = 0; i < idx_lim[0]; i++)
            {
                for (j = 0; j < idx_lim[1]; j++)
                {
                    out[k]+=in[i][j];
                }
            }
        }

//             for (auto& elem : A) {
//                 cnt++;
//                 if (
//                 elem *= multiplier;
// 
//             }
        return out;
    }
        
/*------------------- sum_array_01 END -------------------*/

/*------------------- sum_array_2 -------------------*/
    auto sum_array_2(matlab::data::Array in, int sum_dim)
    {
        using namespace matlab::data;
        ArrayFactory factory;
        std::vector<size_t,2> mat_dim{
                            in1.getDimensions()[0],
                            in1.getDimensions()[1]}
            auto out = factory.createArray<std::complex<double>>(mat_dim);
            auto idx_lim = in.getDimensions();
            for (i = 0; i < idx_lim[0]; i++)
            {
                out[i][j] = 0;
                for (j = 0; j < idx_lim[1]; j++)
                {
                    for (k = 0; k < idx_lim[2]; k++)
                    {
                        out[i][j]+=in[i][j][k];
                    }
                }
            }

//             for (auto& elem : A) {
//                 cnt++;
//                 if (
//                 elem *= multiplier;
// 
//             }
    }
/*------------------- sum_array_2 END-------------------*/
    
    
    
    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        if (inputs.size() != 2) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Two inputs required") }));
        }

        if (inputs[0].getNumberOfElements() != 1) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input multiplier must be a scalar") }));
        }
        
        if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[0].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input multiplier must be a noncomplex scalar double") }));
        }

        if (inputs[1].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input matrix must be type double") }));
        }

        if (inputs[1].getDimensions().size() != 2) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input must be m-by-n dimension") }));
        }
    }
};