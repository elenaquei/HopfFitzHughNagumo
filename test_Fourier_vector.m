% main to test the basic functionality of Fourier_vector

number_of_tests = 0;
Sum = 0;

try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector([1,2,3]);
    Sum = Sum+1;
catch
    warning('a')
end

try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector([1,2,3]');
    if v.nodes ==1
        Sum = Sum+1;
    else
        warning('b')
    end
    number_of_tests = number_of_tests +1;
    w = pad(v,3);
    if w.nodes == 3 && all(w.vector==[0,0,1,2,3,0,0]')
        Sum = Sum+1;
    else
        warning('c')
    end
catch
    warning('d')
end

try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector(rand(1,1,5));
    if all(size(v.vector)==[5,1])
        Sum = Sum+1;
    else
        warning('e')
    end
catch
    warning('f')
end


try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector([1,2,3]);
    w = Fourier_vector([1,2,3,4,5]);
    [w_new,v_new] = padVec(w,v);
    if w_new.nodes == 2 && v_new.nodes ==2
        if all(w_new ==w) && all(v_new.vector == [0,1,2,3,0])
            Sum = Sum+1;
        else
            warning('g')
        end
    end
catch
    warning('h')
end

try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector([1,2,3]);
    w = Fourier_vector([9,0,9]);
    test_vec = Fourier_vector(conv([1,2,3],[9,0,9]));
    vw = prod(v,w);
    if vw == test_vec %|| eq_approx(vw, test_vec)
        Sum = Sum + 1;
    end
end

try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector([1,2,3]);
    w = Fourier_vector([9,0,9]);
    test_vec = Fourier_vector(conv([1,2,3],[9,0,9],'same'));
    vw = prod(v,w,'same');
    if vw == test_vec %|| eq_approx(vw, test_vec)
        Sum = Sum + 1;
    end
end


try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector([1,2,3]);
    w = Fourier_vector([1,2,3]);
    if v==w %|| eq_approx(vw, test_vec)
        Sum = Sum + 1;
    end
end

try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector([1,2,3]);
    w = Fourier_vector([1,2,3]+eps);
    if ~(v==w) && eq_approx(v, w, 10^-8) && eq_approx(v, w)
        Sum = Sum + 1;
    end
end

try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector([1,2,3]);
    w = Fourier_vector([9,0,9]);
    test_vec = Fourier_vector([1,2,3]+[9,0,9]);
    vw = plus(v,w);
    if vw == test_vec
        Sum = Sum + 1;
    end
end


try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector([1,2,3]);
    w = Fourier_vector([0,9,0,9,0]);
    test_vec = Fourier_vector([0,1,2,3,0]+[0,9,0,9,0]);
    vw = plus(v,w);
    if vw == test_vec
        Sum = Sum + 1;
    end
end

try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector(2);
    w = Fourier_vector([0,9,0,9,0]);
    test_vec = Fourier_vector([0,0,2,0,0]+[0,9,0,9,0]);
    vw = plus(v,w);
    if vw == test_vec && vw == plus(w,2)
        Sum = Sum + 1;
    end
end


try 
    number_of_tests = number_of_tests +1;
    v = Fourier_vector(2);
    w = Fourier_vector([0,9,0,9,0]);
    test_vec = Fourier_vector(conv([2],[0,9,0,9,0]));
    vw = prod(v,w);
    if vw == test_vec && vw == prod(w,2)
        Sum = Sum + 1;
    end
end

try 
    number_of_tests = number_of_tests +1;
    w = Fourier_vector([0,-9i,0,9i,0]);
    z = ifft(w);
    if z == ifft([0,9i,0,0,-9i])
        Sum = Sum +1;
    end
end


try 
    number_of_tests = number_of_tests +1;
    w = Fourier_vector([0,9i,0,-9i,0]);
    plot(w,'LineWidth',3)
    Sum = Sum +1;
end


if Sum~=number_of_tests
    error('Somethign is wrong')
else
    fprintf('Success => all tests are good\n')
end