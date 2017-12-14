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
