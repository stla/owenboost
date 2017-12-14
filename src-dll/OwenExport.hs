{-# LANGUAGE ForeignFunctionInterface #-}
module OwenExport
  where
import           Foreign
import           Foreign.C
import           Owen

foreign export ccall owenQexport :: Ptr CInt -> Ptr CDouble -> Ptr CDouble ->
                                 Ptr CDouble -> Ptr CInt -> Ptr CDouble -> IO ()
owenQexport :: Ptr CInt -> Ptr CDouble -> Ptr CDouble -> Ptr CDouble ->
                                                Ptr CInt -> Ptr CDouble -> IO ()
owenQexport nu t delta r n result = do
  nu <- peek nu
  t <- peek t
  n <- peek n
  delta <- peekArray (fromIntegral n) delta
  r <- peekArray (fromIntegral n) r
  (>>=) (owenQ128 nu t delta r) (pokeArray result)

foreign export ccall powenExport :: Ptr CInt -> Ptr CDouble -> Ptr CDouble ->
                  Ptr CDouble -> Ptr CDouble -> Ptr CInt -> Ptr CDouble -> IO ()
powenExport :: Ptr CInt -> Ptr CDouble -> Ptr CDouble -> Ptr CDouble ->
                                Ptr CDouble -> Ptr CInt -> Ptr CDouble -> IO ()
powenExport nu t1 t2 delta1 delta2 n result = do
  nu <- peek nu
  t1 <- peek t1
  t2 <- peek t2
  n <- peek n
  delta1 <- peekArray (fromIntegral n) delta2
  delta2 <- peekArray (fromIntegral n) delta2
  (>>=) (powen nu t1 t2 delta1 delta2) (pokeArray result)
