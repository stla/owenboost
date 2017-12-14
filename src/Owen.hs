{-# LANGUAGE CPP, ForeignFunctionInterface #-}
{-# LANGUAGE QuasiQuotes              #-}
{-# LANGUAGE TemplateHaskell          #-}

module Owen
  where
import           Data.Monoid                  ((<>))
import qualified Data.Vector.Storable         as V
import qualified Data.Vector.Storable.Mutable as VM
import           Foreign
import           Foreign.C.Types
import qualified Language.C.Inline            as C
import qualified Language.C.Inline.Cpp        as CPP

CPP.context (CPP.cppCtx <> C.vecCtx)

CPP.include "<owen128.cpp>"

owenT :: CDouble -> CDouble -> CDouble
owenT h a = [CPP.pure| double {
      owent($(double h), $(double a))
    } |]

studentCDF :: CDouble -> CInt -> [CDouble] -> IO [CDouble]
studentCDF q nu delta = do
  mdelta <- V.thaw (V.fromList delta)
  ptr <- [CPP.exp| double* {
      studentCDF($(double q), $(int nu), $vec-ptr:(double* mdelta), $vec-len:mdelta)
    } |]
  peekArray (length delta) ptr

owenQ :: CInt -> CDouble -> [CDouble] -> [CDouble] -> IO [CDouble]
owenQ nu t delta r = do
  mdelta <- V.thaw (V.fromList delta)
  mr <- V.thaw (V.fromList r)
  ptr <- [CPP.exp| double* {
      owenQ($(int nu), $(double t), $vec-ptr:(double* mdelta), $vec-ptr:(double* mr), $vec-len:mdelta)
    } |]
  VM.clear mdelta
  VM.clear mr
  peekArray (length delta) ptr

owenQ128 :: CInt -> CDouble -> [CDouble] -> [CDouble] -> IO [CDouble]
owenQ128 nu t delta r = do
  mdelta <- V.thaw (V.fromList delta)
  mr <- V.thaw (V.fromList r)
  ptr <- [CPP.exp| double* {
      owenQ128($(int nu), $(double t), $vec-ptr:(double* mdelta), $vec-ptr:(double* mr), $vec-len:mdelta)
    } |]
  VM.clear mdelta
  VM.clear mr
  peekArray (length delta) ptr

powen :: CInt -> CDouble -> CDouble -> [CDouble] -> [CDouble] -> IO [CDouble]
powen nu t1 t2 delta1 delta2 = do
  mdelta1 <- V.thaw (V.fromList delta1)
  mdelta2 <- V.thaw (V.fromList delta2)
  ptr <- [CPP.exp| double* {
      powen4($(int nu), $(double t1), $(double t2), $vec-ptr:(double* mdelta1), $vec-ptr:(double* mdelta2), $vec-len:mdelta1)
    } |]
  VM.clear mdelta1
  VM.clear mdelta2
  peekArray (length delta1) ptr
