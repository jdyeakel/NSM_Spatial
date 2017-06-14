ipbc = function(x,L)
  if x <= (L+2)^2
    check = 0
    if (x <= L+2)
      new_x = x + L*(L+2)
      check = 1
    end
    if (x >= ((L+2)^2 - (L+1)))
      new_x = x - L*(L+2)
      check = 1
    end
    if (x%(L+2) == 1)
      new_x = x + L
      check = 1
    end
    if (x%(L+2) == 0)
      new_x = x - L
      check = 1
    end
    if ((x <= L+2) && (x%(L+2) == 1))
      new_x = x + L*(L+2) + L
      check = 1
    end
    if ((x <= L+2) && (x%(L+2) == 0))
      new_x = x + L*(L+2) - L
      check = 1
    end
    if ((x >= ((L+2)^2 - (L+1))) && (x%(L+2) == 1))
      new_x = x - L*(L+2) + L
      check = 1
    end
    if ((x >= ((L+2)^2 - (L+1))) && (x%(L+2) == 0))
      new_x = x - L*(L+2) - L
      check = 1
    end
    if (check == 0)
      new_x = x
    end
    nn = [new_x+1,new_x-1,new_x+(L+2),new_x-(L+2)]
    return(nn)
  end
end
