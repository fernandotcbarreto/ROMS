### OPERATIONAL WORKING FLOW

countdays=7

while :; do


  current_time=$(date +%H:%M)

  echo $current_time

  if [[ $current_time > "20:00" ]] && [[ $current_time < "20:30" ]]; then 
#   echo 5555

    if [[ $countdays = 7 ]]; then
#   echo 88
      echo $countdays
    
     ./hindcast_forecast.bash
   
     countdays=1
   
    else
#    echo 99

      echo $countdays

      ./forecast_Loop.bash
    
      countdays=$(($countdays + 1))
   
    fi   
   
  fi

sleep 20m


done